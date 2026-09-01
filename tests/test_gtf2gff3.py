"""Regression tests for safe GTF-to-GFF3 conversion."""

import subprocess
import sys
from collections import Counter
from pathlib import Path
from shutil import which
from types import SimpleNamespace

import gff_utils
import pytest
from gtf2gff3 import escape_gff3_component


SCRIPT = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "gtf2gff3.py"
DEFAULT_LOCUS_TAG = object()


def gtf_row(
    feature,
    start,
    stop,
    gene_id,
    strand="+",
    extra="",
    seq_name="chr1",
    score=".",
    phase=None,
    locus_tag=DEFAULT_LOCUS_TAG,
):
    if locus_tag is DEFAULT_LOCUS_TAG:
        locus_tag = f"{gene_id}_tag"
    attributes = f'gene_id "{gene_id}";'
    if locus_tag is not None:
        attributes += f' locus_tag "{locus_tag}";'
    attributes = f"{attributes} {extra}".strip()
    if phase is None:
        phase = "0" if feature.lower() == "cds" else "."
    return "\t".join(
        [
            seq_name,
            "test",
            feature,
            str(start),
            str(stop),
            score,
            strand,
            phase,
            attributes,
        ]
    )


def parse_output(path):
    records = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        assert len(fields) == 9
        attribute_pairs = gff_utils.split_attributes(fields[8])
        records.append(
            {
                "seq_name": fields[0],
                "feature": fields[2],
                "start": int(fields[3]),
                "stop": int(fields[4]),
                "score": fields[5],
                "strand": fields[6],
                "phase": fields[7],
                "raw_attributes": fields[8],
                "attribute_pairs": attribute_pairs,
                "attributes": gff_utils.parse_attributes(fields[8]),
            }
        )
    return records


def conversion_output(tmp_path, lines):
    tmp_path.mkdir(parents=True, exist_ok=True)
    annotation = tmp_path / "input.gtf"
    annotation.write_text("\n".join(lines) + "\n")
    output = tmp_path / "nested" / "output.gff"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-a", str(annotation), "-o", str(output)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert output.is_file() and output.stat().st_size > 0
    assert output.read_text().startswith("##gff-version 3\n")
    return output


def run_conversion(tmp_path, lines):
    return parse_output(conversion_output(tmp_path, lines))


def test_gene_rna_and_cds_keep_their_own_coordinates(tmp_path):
    records = run_conversion(
        tmp_path,
        [
            gtf_row("gene", 100, 500, "shared", phase="1"),
            gtf_row("CDS", 130, 420, "shared", phase="2"),
            gtf_row("rRNA", 180, 260, "shared", phase="0"),
        ],
    )

    by_feature = {record["feature"]: record for record in records}
    assert (by_feature["gene"]["start"], by_feature["gene"]["stop"]) == (100, 500)
    assert (by_feature["CDS"]["start"], by_feature["CDS"]["stop"]) == (130, 420)
    assert (by_feature["rRNA"]["start"], by_feature["rRNA"]["stop"]) == (180, 260)
    assert by_feature["gene"]["phase"] == "."
    assert by_feature["CDS"]["phase"] == "2"
    assert by_feature["rRNA"]["phase"] == "."
    assert (
        by_feature["CDS"]["attributes"]["parent"]
        == by_feature["gene"]["attributes"]["id"]
    )
    assert (
        by_feature["rRNA"]["attributes"]["parent"]
        == by_feature["gene"]["attributes"]["id"]
    )


def test_rna_without_gene_gets_a_parent_and_does_not_crash(tmp_path):
    records = run_conversion(tmp_path, [gtf_row("tRNA", 200, 275, "rna_only")])
    assert [record["feature"] for record in records] == ["gene", "tRNA"]
    assert [(record["start"], record["stop"]) for record in records] == [
        (200, 275),
        (200, 275),
    ]
    assert records[1]["attributes"]["parent"] == records[0]["attributes"]["id"]


def test_orphan_cds_and_rna_share_one_synthetic_parent(tmp_path):
    records = run_conversion(
        tmp_path,
        [
            gtf_row("CDS", 100, 150, "mixed_orphan"),
            gtf_row("rRNA", 120, 140, "mixed_orphan"),
        ],
    )

    by_feature = {record["feature"]: record for record in records}
    assert set(by_feature) == {"gene", "CDS", "rRNA"}
    assert (by_feature["gene"]["start"], by_feature["gene"]["stop"]) == (100, 150)
    assert (by_feature["CDS"]["start"], by_feature["CDS"]["stop"]) == (100, 150)
    assert (by_feature["rRNA"]["start"], by_feature["rRNA"]["stop"]) == (120, 140)
    assert (
        by_feature["CDS"]["attributes"]["parent"]
        == by_feature["gene"]["attributes"]["id"]
    )
    assert (
        by_feature["rRNA"]["attributes"]["parent"]
        == by_feature["gene"]["attributes"]["id"]
    )


@pytest.mark.parametrize("strand", ["+", "-"])
def test_orphan_cds_parent_matches_child_bounds_and_preserves_columns(tmp_path, strand):
    records = run_conversion(
        tmp_path,
        [
            gtf_row(
                "CDS",
                100,
                150,
                "cds_only",
                strand=strand,
                score="5",
                phase="2",
            )
        ],
    )
    by_feature = {record["feature"]: record for record in records}
    assert (by_feature["gene"]["start"], by_feature["gene"]["stop"]) == (100, 150)
    assert (by_feature["CDS"]["start"], by_feature["CDS"]["stop"]) == (100, 150)
    assert by_feature["gene"]["strand"] == strand
    assert by_feature["CDS"]["strand"] == strand
    assert (by_feature["gene"]["score"], by_feature["gene"]["phase"]) == (".", ".")
    assert (by_feature["CDS"]["score"], by_feature["CDS"]["phase"]) == ("5", "2")


def test_multiple_orphan_cds_records_are_all_preserved(tmp_path):
    records = run_conversion(
        tmp_path,
        [
            gtf_row("CDS", 100, 150, "split_cds"),
            gtf_row("CDS", 200, 249, "split_cds", phase="1"),
        ],
    )
    cds = [record for record in records if record["feature"] == "CDS"]
    gene = next(record for record in records if record["feature"] == "gene")
    assert [(record["start"], record["stop"]) for record in cds] == [
        (100, 150),
        (200, 249),
    ]
    assert (gene["start"], gene["stop"]) == (100, 249)
    assert [record["phase"] for record in cds] == ["0", "1"]


def test_unknown_only_gtf_is_preserved_with_its_own_identifier(tmp_path):
    records = run_conversion(
        tmp_path,
        [gtf_row("repeat_region", 700, 730, "repeat_one", phase="2")],
    )
    assert len(records) == 1
    assert records[0]["feature"] == "repeat_region"
    assert records[0]["attributes"]["id"] == "misc1"
    assert records[0]["phase"] == "."
    assert records[0]["attributes"]["locus_tag"] == "repeat_one_tag"
    assert records[0]["attributes"]["gene_id"] == "repeat_one"


def test_explicit_locus_tag_wins_and_missing_tag_falls_back_to_gene_id(tmp_path):
    explicit = run_conversion(
        tmp_path / "explicit",
        [
            gtf_row("gene", 100, 200, "source_gene", locus_tag="source_locus"),
            gtf_row("CDS", 120, 180, "source_gene", locus_tag="source_locus"),
        ],
    )
    fallback = run_conversion(
        tmp_path / "fallback",
        [
            gtf_row("gene", 300, 400, "fallback_gene", locus_tag=None),
            gtf_row("CDS", 320, 380, "fallback_gene", locus_tag=""),
        ],
    )

    assert {record["attributes"]["locus_tag"] for record in explicit} == {
        "source_locus"
    }
    assert {record["attributes"]["gene_id"] for record in explicit} == {"source_gene"}
    assert {record["attributes"]["locus_tag"] for record in fallback} == {
        "fallback_gene"
    }
    assert {record["attributes"]["gene_id"] for record in fallback} == {"fallback_gene"}


def test_reannotation_downstream_receives_source_gene_and_locus_tag(tmp_path):
    pytest.importorskip("pandas")
    from reannotate_orfs import generate_annotation_dict

    output = conversion_output(
        tmp_path,
        [
            gtf_row("gene", 100, 200, "source_gene", locus_tag="source_locus"),
            gtf_row("CDS", 120, 180, "source_gene", locus_tag="source_locus"),
        ],
    )

    downstream = generate_annotation_dict(SimpleNamespace(annotation_path=output))
    gene_id, locus_tag, *_ = downstream["chr1:120-180:+"]
    assert gene_id == "source_gene"
    assert locus_tag == "source_locus"


def test_generated_attributes_are_unique_and_values_are_gff3_encoded(tmp_path):
    output = conversion_output(
        tmp_path,
        [
            gtf_row(
                "gene",
                100,
                500,
                "canonical tag",
                extra=(
                    'ID "stale_gene"; Parent "stale_parent"; Gene_ID "stale"; '
                    'Note "a=b;c&d,e already%3Dencoded bare%percent / résumé"; '
                    'custom=key "résumé / draft";'
                ),
            ),
            gtf_row(
                "CDS",
                130,
                420,
                "canonical tag",
                extra='ID "stale_cds"; Parent "stale_parent";',
            ),
        ],
    )
    records = parse_output(output)
    gene = next(record for record in records if record["feature"] == "gene")
    cds = next(record for record in records if record["feature"] == "CDS")

    for record in records:
        key_counts = Counter(key.lower() for key, _ in record["attribute_pairs"])
        assert key_counts["id"] == 1
        assert key_counts["locus_tag"] == 1
        assert key_counts["gene_id"] == 1
    assert Counter(key.lower() for key, _ in gene["attribute_pairs"])["parent"] == 0
    assert Counter(key.lower() for key, _ in cds["attribute_pairs"])["parent"] == 1

    assert gene["attributes"]["locus_tag"] == "canonical tag_tag"
    assert gene["attributes"]["gene_id"] == "canonical tag"
    assert gene["attributes"]["note"] == (
        "a%3Db%3Bc%26d%2Ce already%253Dencoded bare%25percent / résumé"
    )
    assert gene["attributes"]["custom%3dkey"] == "résumé / draft"
    assert cds["attributes"]["parent"] == gene["attributes"]["id"]

    genome_tools = which("gt")
    if genome_tools is not None:
        result = subprocess.run(
            [genome_tools, "gff3validator", str(output)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("a=b;c&d,e%", "a%3Db%3Bc%26d%2Ce%25"),
        ("\x00\t\n\r\x1f\x7f", "%00%09%0A%0D%1F%7F"),
        ('café / alpha "beta"', 'café / alpha "beta"'),
        ("%3D %3d %00 %C3%A9", "%253D %253d %2500 %25C3%25A9"),
        ("% %2 %GG %%3B", "%25 %252 %25GG %25%253B"),
    ],
)
def test_gtf_component_escaping_encodes_every_literal_percent(value, expected):
    assert escape_gff3_component(value) == expected


def test_gff3_passthrough_keeps_existing_percent_escapes(tmp_path):
    lines = [
        "##gff-version 3",
        "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1;Note=already%3Dencoded;",
    ]
    output = conversion_output(tmp_path, lines)
    assert output.read_text() == "\n".join(lines) + "\n"


def test_shuffled_gtf_rows_produce_identical_output_and_ids(tmp_path):
    lines = [
        gtf_row("CDS", 300, 360, "alpha", phase="2"),
        gtf_row("rRNA", 220, 250, "alpha", phase="1"),
        gtf_row("gene", 100, 400, "alpha", phase="2"),
        gtf_row("CDS", 120, 180, "alpha", phase="1"),
        gtf_row("rRNA", 600, 630, "beta", phase="2"),
        gtf_row("CDS", 500, 650, "beta", phase="0"),
        gtf_row("repeat_region", 900, 930, "repeat_b", phase="1"),
        gtf_row("repeat_region", 800, 830, "repeat_a", phase="2"),
    ]
    shuffled = [lines[index] for index in (6, 3, 5, 1, 7, 0, 4, 2)]

    first = conversion_output(tmp_path / "first", lines)
    second = conversion_output(tmp_path / "second", shuffled)
    assert first.read_text() == second.read_text()

    records = parse_output(first)
    cds_ids = {
        (record["attributes"]["gene_id"], record["start"]): record["attributes"]["id"]
        for record in records
        if record["feature"] == "CDS"
    }
    assert cds_ids == {
        ("alpha", 120): "cds1",
        ("alpha", 300): "cds2",
        ("beta", 500): "cds3",
    }
    misc_ids = [
        record["attributes"]["id"]
        for record in records
        if record["feature"] == "repeat_region"
    ]
    assert misc_ids == ["misc1", "misc2"]


def run_invalid(tmp_path, lines):
    annotation = tmp_path / "invalid.gtf"
    output = tmp_path / "must-not-exist.gff"
    annotation.write_text("\n".join(lines) + "\n")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-a", str(annotation), "-o", str(output)],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert not output.exists()
    assert "Traceback" not in result.stderr
    return result.stderr


def test_conflicting_locus_tags_for_one_gene_id_are_rejected(tmp_path):
    stderr = run_invalid(
        tmp_path,
        [
            gtf_row("gene", 100, 200, "conflicting", locus_tag="first"),
            gtf_row("CDS", 120, 180, "conflicting", locus_tag="second"),
        ],
    )
    assert "conflicting locus_tag values" in stderr
    assert "'first'" in stderr
    assert "'second'" in stderr


def test_duplicate_gene_records_are_rejected_without_replacing_existing_output(
    tmp_path,
):
    annotation = tmp_path / "duplicate.gtf"
    annotation.write_text(
        "\n".join(
            [
                gtf_row("gene", 100, 200, "duplicate"),
                gtf_row("gene", 100, 200, "duplicate"),
            ]
        )
        + "\n"
    )
    output = tmp_path / "nested" / "existing.gff"
    output.parent.mkdir()
    output.write_text("existing output\n")

    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-a", str(annotation), "-o", str(output)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "duplicate gene records" in result.stderr
    assert "Traceback" not in result.stderr
    assert output.read_text() == "existing output\n"
    assert not list(output.parent.glob(f".{output.name}.*.tmp"))


@pytest.mark.parametrize(
    ("second_row", "message"),
    [
        (
            gtf_row("CDS", 130, 160, "conflict", seq_name="chr2"),
            "spans multiple sequence IDs",
        ),
        (
            gtf_row("CDS", 130, 160, "conflict", strand="-"),
            "spans multiple strands",
        ),
    ],
)
def test_same_gene_id_cannot_span_sequence_or_strand(tmp_path, second_row, message):
    stderr = run_invalid(
        tmp_path,
        [gtf_row("CDS", 100, 120, "conflict"), second_row],
    )
    assert message in stderr


def test_child_outside_explicit_gene_is_rejected(tmp_path):
    stderr = run_invalid(
        tmp_path,
        [
            gtf_row("gene", 100, 200, "bad_bounds"),
            gtf_row("CDS", 90, 150, "bad_bounds"),
        ],
    )
    assert "lies outside gene_id 'bad_bounds'" in stderr


def test_mixed_gff3_and_gtf_fails_with_a_clear_message(tmp_path):
    stderr = run_invalid(
        tmp_path,
        [
            gtf_row("gene", 100, 200, "gtf_gene"),
            "chr1\ttest\tgene\t300\t400\t.\t+\t.\tID=gff_gene;locus_tag=gff_gene;",
        ],
    )
    assert "mixed GFF3 and GTF attribute syntax" in stderr


def test_unsupported_identifier_style_fails_instead_of_omitting_output(tmp_path):
    stderr = run_invalid(
        tmp_path,
        ['chr1\ttest\tgene\t100\t200\t.\t+\t.\tgene_name "nameless";'],
    )
    assert "GTF input contains no gene_id attribute" in stderr
