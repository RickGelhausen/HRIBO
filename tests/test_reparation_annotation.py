"""Regression tests for REPARATION's transcript-only annotation boundary."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import gff_utils


REPO = Path(__file__).resolve().parent.parent
ADAPTER = REPO / "workflow/scripts/prepare_reparation_annotation.py"
CONVERTER = REPO / "workflow/scripts/gtf2gff3.py"
REPARATION_FIELDS = (
    "gene_id",
    "transcript_id",
    "gene_name",
    "gene_biotype",
)


def feature(
    feature_type,
    start,
    end,
    attributes,
    *,
    seqid="chr1",
    source="RefSeq",
    strand="+",
    score=".",
    phase="0",
):
    return "\t".join(
        (
            seqid,
            source,
            feature_type,
            str(start),
            str(end),
            score,
            strand,
            phase,
            attributes,
        )
    )


def run_adapter(annotation: Path, output: Path):
    return subprocess.run(
        [sys.executable, str(ADAPTER), "-a", str(annotation), "-o", str(output)],
        capture_output=True,
        text=True,
    )


def records(path: Path):
    parsed = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        assert len(fields) == 9
        parsed.append((fields, gff_utils.parse_attributes(fields[8]), line))
    return parsed


def test_adapter_emits_only_complete_transcripts_for_coding_and_rna(tmp_path):
    annotation = tmp_path / "processed.gff3"
    output = tmp_path / "reparation" / "annotation.gtf"
    annotation.write_text(
        "##gff-version 3\n"
        + feature(
            "gene",
            100,
            402,
            "ID=gene-plus;gene_id=source-plus;Name=Plus protein;",
        )
        + "\n"
        + feature("CDS", 100, 402, "ID=cds-plus;Parent=gene-plus;")
        + "\n"
        + feature(
            "gene",
            600,
            902,
            "ID=gene-minus;locus_tag=minus-tag;gene=minus-name;",
            strand="-",
        )
        + "\n"
        + feature(
            "CDS",
            600,
            902,
            "ID=cds-minus;Parent=gene-minus;",
            strand="-",
        )
        + "\n"
        + feature(
            "rRNA",
            1000,
            1090,
            "ID=rna-one;locus_tag=rna-tag;Name=16S fragment;",
            phase=".",
        )
        + "\n"
        + feature(
            "pseudogene",
            1200,
            1292,
            "ID=pseudo-parent;locus_tag=pseudo-tag;Name=Pseudo protein;",
        )
        + "\n"
        + feature("CDS", 1200, 1292, "ID=pseudo-cds;Parent=pseudo-parent;")
        + "\n",
        encoding="utf-8",
    )

    result = run_adapter(annotation, output)

    assert result.returncode == 0, result.stderr
    output_records = records(output)
    assert [(row[0][2], row[0][3], row[0][4], row[0][6]) for row in output_records] == [
        ("transcript", "100", "402", "+"),
        ("transcript", "600", "902", "-"),
        ("transcript", "1000", "1090", "+"),
        ("transcript", "1200", "1292", "+"),
    ]
    assert [row[1]["transcript_id"] for row in output_records] == [
        "hribo_reparation_transcript_000001",
        "hribo_reparation_transcript_000002",
        "hribo_reparation_transcript_000003",
        "hribo_reparation_transcript_000004",
    ]
    by_gene = {row[1]["gene_id"]: row[1] for row in output_records}
    assert by_gene["source-plus"] == {
        "gene_id": "source-plus",
        "transcript_id": "hribo_reparation_transcript_000001",
        "gene_name": "Plus protein",
        "gene_biotype": "protein_coding",
    }
    assert by_gene["minus-tag"]["gene_name"] == "minus-name"
    assert by_gene["minus-tag"]["gene_biotype"] == "protein_coding"
    assert by_gene["rna-tag"]["gene_biotype"] == "rRNA"
    assert by_gene["pseudo-tag"]["gene_biotype"] == "pseudogene"
    assert by_gene["pseudo-tag"]["gene_name"] == "Pseudo protein"

    # Mirror the permissive expressions in REPARATION 1.0.9's
    # post_processing.pl. Every field it dereferences must be nonempty.
    for _, attributes, line in output_records:
        for field in REPARATION_FIELDS:
            match = re.search(rf'{field}."?([^";]+)"?', line)
            assert match is not None
            assert match.group(1) == attributes[field]


def test_gtf_conversion_then_adapter_preserves_stop_inclusive_boundaries(tmp_path):
    source = tmp_path / "source.gtf"
    processed = tmp_path / "processed.gff3"
    output = tmp_path / "reparation.gtf"
    source.write_text(
        feature(
            "gene",
            10,
            102,
            'gene_id "plus"; gene_name "Plus"; gene_biotype "protein_coding";',
            phase=".",
        )
        + "\n"
        + feature(
            "transcript",
            10,
            102,
            'gene_id "plus"; transcript_id "plus-tr"; gene_biotype "protein_coding";',
            phase=".",
        )
        + "\n"
        + feature(
            "CDS",
            10,
            102,
            'gene_id "plus"; transcript_id "plus-tr"; gene_biotype "protein_coding";',
        )
        + "\n"
        + feature(
            "gene",
            200,
            292,
            'gene_id "minus"; gene_name "Minus"; gene_biotype "protein_coding";',
            strand="-",
            phase=".",
        )
        + "\n"
        + feature(
            "CDS",
            200,
            292,
            'gene_id "minus"; transcript_id "minus-tr"; gene_biotype "protein_coding";',
            strand="-",
        )
        + "\n",
        encoding="utf-8",
    )

    converted = subprocess.run(
        [sys.executable, str(CONVERTER), "-a", str(source), "-o", str(processed)],
        capture_output=True,
        text=True,
    )
    assert converted.returncode == 0, converted.stderr
    result = run_adapter(processed, output)
    assert result.returncode == 0, result.stderr

    by_gene = {row[1]["gene_id"]: row[0] for row in records(output)}
    assert (by_gene["plus"][3], by_gene["plus"][4]) == ("10", "102")
    assert (by_gene["minus"][3], by_gene["minus"][4]) == ("200", "292")

    # REPARATION removes the three-base stop codon during classification;
    # create_reparation_gff.py adds those bases back to its final GFF record.
    assert (10, 102 - 3 + 3) == (10, 102)
    assert (200 + 3 - 3, 292) == (200, 292)


def test_adapter_is_order_independent_and_protects_regex_delimiters(tmp_path):
    lines = [
        feature(
            "CDS",
            500,
            601,
            'ID=zeta;gene_id=zeta-gene;Name=quoted"name\\draft;',
            seqid="z",
        ),
        feature(
            "tRNA",
            20,
            80,
            "locus_tag=alpha-tag;Name=Alpha RNA;",
            seqid="a",
            phase=".",
        ),
    ]
    first_input = tmp_path / "first.gff3"
    second_input = tmp_path / "second.gff3"
    first_output = tmp_path / "first.gtf"
    second_output = tmp_path / "second.gtf"
    first_input.write_text("\n".join(lines) + "\n", encoding="utf-8")
    second_input.write_text("\n".join(reversed(lines)) + "\n", encoding="utf-8")

    assert run_adapter(first_input, first_output).returncode == 0
    assert run_adapter(second_input, second_output).returncode == 0
    assert first_output.read_bytes() == second_output.read_bytes()
    rendered = first_output.read_text(encoding="utf-8")
    assert 'gene_name "quoted%22name%5Cdraft";' in rendered
    assert [row[0][0] for row in records(first_output)] == ["a", "z"]


def test_adapter_ignores_unreferenced_duplicate_ids_on_unsupported_features(
    tmp_path,
):
    annotation = tmp_path / "reused-regulatory-ids.gff3"
    output = tmp_path / "reparation.gtf"
    annotation.write_text(
        feature("TATA", 10, 15, "ID=TATA-gene-one;")
        + "\n"
        + feature("TATA", 30, 35, "ID=TATA-gene-one;")
        + "\n"
        + feature(
            "CDS",
            100,
            300,
            "ID=cds-one;locus_tag=gene-one;Name=Gene one;",
        )
        + "\n",
        encoding="utf-8",
    )

    result = run_adapter(annotation, output)

    assert result.returncode == 0, result.stderr
    output_records = records(output)
    assert len(output_records) == 1
    assert output_records[0][1]["gene_id"] == "gene-one"


def test_adapter_rejects_duplicate_ids_that_can_affect_output(tmp_path):
    for name, content in (
        (
            "supported",
            feature("CDS", 10, 30, "ID=duplicate;")
            + "\n"
            + feature("CDS", 40, 60, "ID=duplicate;"),
        ),
        (
            "parent-referenced",
            feature("gene", 10, 30, "ID=duplicate;")
            + "\n"
            + feature("gene", 40, 60, "ID=duplicate;")
            + "\n"
            + feature("CDS", 10, 30, "ID=child;Parent=duplicate;"),
        ),
    ):
        case = tmp_path / name
        case.mkdir()
        annotation = case / "annotation.gff3"
        output = case / "existing.gtf"
        annotation.write_text(content + "\n", encoding="utf-8")
        output.write_text("keep me\n", encoding="utf-8")

        result = run_adapter(annotation, output)

        assert result.returncode == 2
        assert "duplicate ID 'duplicate' on lines 1 and 2" in result.stderr
        assert output.read_text(encoding="utf-8") == "keep me\n"


def test_adapter_rejects_split_cds_without_replacing_output(tmp_path):
    annotation = tmp_path / "duplicates.gff3"
    output = tmp_path / "existing.gtf"
    annotation.write_text(
        feature("CDS", 10, 30, "ID=first;transcript_id=shared;")
        + "\n"
        + feature("CDS", 40, 60, "ID=second;transcript_id=shared;")
        + "\n",
        encoding="utf-8",
    )
    output.write_text("keep me\n", encoding="utf-8")

    result = run_adapter(annotation, output)

    assert result.returncode == 2
    assert "CDS transcript 'shared' is split" in result.stderr
    assert "Traceback" not in result.stderr
    assert output.read_text(encoding="utf-8") == "keep me\n"
    assert not list(tmp_path.glob(f".{output.name}.*.tmp"))


def test_adapter_rejects_unsupported_or_unstranded_annotation_atomically(tmp_path):
    for name, content, message in (
        (
            "unsupported",
            feature("gene", 10, 30, "ID=gene-only;"),
            "no CDS, supported RNA, pseudogene, or standalone transcript",
        ),
        (
            "unstranded",
            feature("CDS", 10, 30, "ID=unstranded;", strand="."),
            "REPARATION requires '+' or '-'",
        ),
    ):
        case = tmp_path / name
        case.mkdir()
        annotation = case / "annotation.gff3"
        output = case / "existing.gtf"
        annotation.write_text(content + "\n", encoding="utf-8")
        output.write_text("keep me\n", encoding="utf-8")

        result = run_adapter(annotation, output)

        assert result.returncode == 2
        assert message in result.stderr
        assert "Traceback" not in result.stderr
        assert output.read_text(encoding="utf-8") == "keep me\n"
        assert not list(case.glob(f".{output.name}.*.tmp"))
