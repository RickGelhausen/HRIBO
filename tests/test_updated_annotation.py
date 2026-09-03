"""Direct regression tests for the deterministic updated-annotation builder."""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "workflow" / "scripts" / "build_updated_annotation.py"
HEADER = "##gff-version 3\n"


def feature(
    source,
    feature_type,
    start,
    end,
    attributes,
    *,
    strand="+",
    phase=".",
):
    return (
        "\t".join(
            [
                "chr",
                source,
                feature_type,
                str(start),
                str(end),
                ".",
                strand,
                phase,
                attributes,
            ]
        )
        + "\n"
    )


def run_builder(annotation, predictions, output):
    command = [sys.executable, str(SCRIPT), "--annotation", str(annotation)]
    if predictions is not None:
        command.extend(["--predictions", *map(str, predictions)])
    command.extend(["--output", str(output)])
    return subprocess.run(command, capture_output=True, text=True, cwd=SCRIPT.parent)


def records(path):
    return [
        line.split("\t")
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def attributes(record):
    return {
        key: value
        for field in record[8].split(";")
        if field
        for key, value in [field.split("=", 1)]
    }


def test_namespaces_colliding_ids_and_rewrites_only_internal_references(tmp_path):
    annotation = tmp_path / "annotation.gff3"
    alpha = tmp_path / "alpha.gff3"
    beta = tmp_path / "beta.gff3"
    output = tmp_path / "updated.gff3"
    reversed_output = tmp_path / "updated-reversed.gff3"

    original_attributes = "Note=keep%3Bexactly;ID=shared;custom=Value=With=Equals;empty="
    annotation.write_text(
        HEADER
        + "# this input comment is intentionally not another version header\n"
        + feature("reference", "gene", 1, 100, original_attributes),
        encoding="utf-8",
    )
    alpha.write_text(
        HEADER
        + feature("alpha", "gene", 200, 280, "ID=shared;Name=alpha;")
        + feature(
            "alpha",
            "mRNA",
            200,
            280,
            "ID=transcript;Parent=shared,external_parent;"
            "Derives_from=shared,beta_only,external_derivation;",
        ),
        encoding="utf-8",
    )
    beta.write_text(
        HEADER
        + feature("beta", "gene", 300, 380, "ID=shared;Name=beta;")
        + feature("beta", "gene", 400, 480, "ID=beta_only;"),
        encoding="utf-8",
    )

    result = run_builder(annotation, [alpha, beta], output)
    assert result.returncode == 0, result.stderr
    reverse_result = run_builder(annotation, [beta, alpha], reversed_output)
    assert reverse_result.returncode == 0, reverse_result.stderr

    assert output.read_bytes() == reversed_output.read_bytes()
    assert output.read_text(encoding="utf-8").count(HEADER) == 1
    output_records = records(output)
    assert len(output_records) == 5

    original = next(record for record in output_records if record[1] == "reference")
    assert original[8] == original_attributes

    by_source_and_type = {
        (record[1], record[2]): attributes(record) for record in output_records
    }
    assert by_source_and_type[("alpha", "gene")]["ID"] == "alpha:shared"
    beta_ids = {
        attributes(record)["ID"]
        for record in output_records
        if record[1] == "beta"
    }
    assert beta_ids == {
        "beta:shared",
        "beta:beta_only",
    }

    transcript = by_source_and_type[("alpha", "mRNA")]
    assert transcript["ID"] == "alpha:transcript"
    assert transcript["Parent"] == "alpha:shared,external_parent"
    assert (
        transcript["Derives_from"]
        == "alpha:shared,beta:beta_only,external_derivation"
    )


@pytest.mark.parametrize("prediction_content", ["", HEADER], ids=["empty", "header-only"])
def test_empty_prediction_file_produces_original_only(tmp_path, prediction_content):
    annotation = tmp_path / "annotation.gff3"
    prediction = tmp_path / "prediction.gff3"
    output = tmp_path / "updated.gff3"
    original_record = feature("reference", "gene", 1, 10, "ID=original;Name=Exact;")
    annotation.write_text(HEADER + original_record, encoding="utf-8")
    prediction.write_text(prediction_content, encoding="utf-8")

    result = run_builder(annotation, [prediction], output)

    assert result.returncode == 0, result.stderr
    assert output.read_text(encoding="utf-8") == HEADER + original_record


def test_predictions_option_may_be_omitted(tmp_path):
    annotation = tmp_path / "annotation.gff3"
    output = tmp_path / "updated.gff3"
    original_record = feature("reference", "gene", 1, 10, "ID=original;")
    annotation.write_text(HEADER + original_record, encoding="utf-8")

    result = run_builder(annotation, None, output)

    assert result.returncode == 0, result.stderr
    assert output.read_text(encoding="utf-8") == HEADER + original_record


def test_existing_gff3_escapes_are_preserved_and_stray_percent_is_encoded(tmp_path):
    annotation = tmp_path / "annotation.gff3"
    prediction = tmp_path / "prediction.gff3"
    output = tmp_path / "updated.gff3"
    annotation.write_text(
        HEADER + feature("reference", "gene", 1, 10, "ID=original;"),
        encoding="utf-8",
    )
    prediction.write_text(
        HEADER
        + feature("caller%2Fv1", "gene", 20, 40, "ID=gene%3A1;")
        + feature(
            "caller%2Fv1",
            "mRNA",
            20,
            40,
            "ID=transcript%2F1;Parent=gene%3A1;",
        )
        + feature("caller%XZ", "gene", 50, 60, "ID=loose%name;"),
        encoding="utf-8",
    )

    result = run_builder(annotation, [prediction], output)

    assert result.returncode == 0, result.stderr
    prediction_attributes = {
        record[2]: attributes(record)
        for record in records(output)
        if record[1] == "caller%2Fv1"
    }
    assert prediction_attributes["gene"]["ID"] == "caller%2Fv1:gene%3A1"
    assert prediction_attributes["mRNA"]["ID"] == "caller%2Fv1:transcript%2F1"
    assert prediction_attributes["mRNA"]["Parent"] == "caller%2Fv1:gene%3A1"
    stray = next(record for record in records(output) if record[1] == "caller%XZ")
    assert attributes(stray)["ID"] == "caller%25XZ:loose%25name"


def test_generated_ids_avoid_reference_and_encoded_namespace_collisions(tmp_path):
    annotation = tmp_path / "annotation.gff3"
    prediction = tmp_path / "prediction.gff3"
    output = tmp_path / "updated.gff3"
    annotation.write_text(
        HEADER
        + feature("reference", "gene", 1, 10, "ID=alpha:shared;")
        + feature("reference", "gene", 20, 30, "ID=caller%2Fv1:gene;"),
        encoding="utf-8",
    )
    prediction.write_text(
        HEADER
        + feature("alpha", "gene", 40, 50, "ID=shared;")
        + feature("alpha", "mRNA", 40, 50, "ID=transcript;Parent=shared;")
        + feature("caller%2Fv1", "gene", 60, 70, "ID=gene;")
        + feature("caller/v1", "gene", 80, 90, "ID=gene;"),
        encoding="utf-8",
    )

    result = run_builder(annotation, [prediction], output)

    assert result.returncode == 0, result.stderr
    output_records = records(output)
    alpha_gene = next(
        record
        for record in output_records
        if record[1] == "alpha" and record[2] == "gene"
    )
    alpha_transcript = next(
        record for record in output_records if record[1] == "alpha" and record[2] == "mRNA"
    )
    assert attributes(alpha_gene)["ID"] == "alpha:shared:prediction"
    assert attributes(alpha_transcript)["Parent"] == "alpha:shared:prediction"

    encoded_collision_ids = {
        record[1]: attributes(record)["ID"]
        for record in output_records
        if record[1] in {"caller%2Fv1", "caller/v1"}
    }
    assert encoded_collision_ids == {
        "caller%2Fv1": "caller%2Fv1:gene:prediction",
        "caller/v1": "caller%2Fv1:gene:prediction2",
    }


def test_malformed_input_reports_context_and_preserves_existing_output(tmp_path):
    annotation = tmp_path / "annotation.gff3"
    malformed = tmp_path / "malformed.gff3"
    output = tmp_path / "updated.gff3"
    stale = HEADER + feature("old", "gene", 1, 10, "ID=stale;")
    annotation.write_text(
        HEADER + feature("reference", "gene", 1, 10, "ID=original;"),
        encoding="utf-8",
    )
    malformed.write_text("chr\tpredictor\tCDS\t1\t10\t.\t+\t0\n", encoding="utf-8")
    output.write_text(stale, encoding="utf-8")

    result = run_builder(annotation, [malformed], output)

    assert result.returncode != 0
    assert f"{malformed}:1" in result.stderr
    assert "expected 9 tab-separated GFF3 columns" in result.stderr
    assert output.read_text(encoding="utf-8") == stale


def test_combined_output_passes_strict_genometools_validation(tmp_path):
    genome_tools = shutil.which("gt")
    assert genome_tools is not None, "GenomeTools is required for strict GFF3 tests"

    annotation = tmp_path / "annotation.gff3"
    alpha = tmp_path / "alpha.gff3"
    beta = tmp_path / "beta.gff3"
    output = tmp_path / "updated.gff3"
    annotation.write_text(
        HEADER
        + feature("reference", "gene", 1, 100, "ID=reference_gene;Name=reference;")
        + feature(
            "reference",
            "mRNA",
            1,
            100,
            "ID=reference_transcript;Parent=reference_gene;",
        )
        + feature(
            "reference",
            "CDS",
            10,
            90,
            "ID=reference_cds;Parent=reference_transcript;",
            phase="0",
        ),
        encoding="utf-8",
    )
    alpha.write_text(
        HEADER
        + feature("alpha", "gene", 200, 300, "ID=prediction_gene;")
        + feature(
            "alpha", "mRNA", 200, 300, "ID=prediction_tx;Parent=prediction_gene;"
        )
        + feature(
            "alpha",
            "CDS",
            210,
            290,
            "ID=prediction_cds;Parent=prediction_tx;",
            phase="0",
        ),
        encoding="utf-8",
    )
    beta.write_text(
        HEADER
        + feature("beta", "gene", 200, 300, "ID=prediction_gene;")
        + feature(
            "beta", "mRNA", 200, 300, "ID=prediction_tx;Parent=prediction_gene;"
        )
        + feature(
            "beta",
            "CDS",
            210,
            290,
            "ID=prediction_cds;Parent=prediction_tx;",
            phase="0",
        ),
        encoding="utf-8",
    )

    result = run_builder(annotation, [alpha, beta], output)
    assert result.returncode == 0, result.stderr
    validation = subprocess.run(
        [genome_tools, "gff3validator", str(output)], capture_output=True, text=True
    )

    assert validation.returncode == 0, validation.stderr
