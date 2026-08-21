"""Tests for the input preflight validation.

Each test asserts on the check identifier rather than on message text, so that
wording can be improved without breaking the suite.
"""

from pathlib import Path

import pytest

from lib import checks
from lib.validation import Severity, parse_annotation, parse_fasta


def check_ids(report, severity=None):
    return {
        f.check for f in report.findings if severity is None or f.severity is severity
    }


# --------------------------------------------------------------------------
# The happy path
# --------------------------------------------------------------------------


def test_valid_inputs_produce_no_errors(genome_file, annotation_file, config, samples):
    genome_report, genome_records = checks.check_genome(genome_file)
    annotation_report, gff_records = checks.check_annotation(annotation_file)
    consistency = checks.check_reference_consistency(genome_records, gff_records)

    for report in (genome_report, annotation_report, consistency):
        assert report.ok, report.to_text()

    assert checks.check_sample_sheet(samples).ok
    assert checks.check_fastq_files(samples).ok
    assert checks.check_config_semantics(config, samples).ok


# --------------------------------------------------------------------------
# Genome / annotation identifier consistency
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "replacement, expected_reason",
    [
        ("NC_000913", "version suffix"),
        ("nc_000913.3", "capitalisation"),
        ("Escherichia", "FASTA description"),
    ],
)
def test_seqid_mismatch_is_detected_and_explained(
    genome_file, annotation_file, replacement, expected_reason
):
    annotation_file.write_text(
        annotation_file.read_text().replace("NC_000913.3\t", f"{replacement}\t")
    )

    _, genome_records = checks.check_genome(genome_file)
    _, gff_records = checks.check_annotation(annotation_file)
    report = checks.check_reference_consistency(genome_records, gff_records)

    assert "REFERENCE_ID_MISMATCH" in check_ids(report, Severity.ERROR)
    finding = next(f for f in report.findings if f.check == "REFERENCE_ID_MISMATCH")
    assert expected_reason in finding.hint, finding.hint


def test_matching_identifiers_are_accepted(genome_file, annotation_file):
    _, genome_records = checks.check_genome(genome_file)
    _, gff_records = checks.check_annotation(annotation_file)
    assert checks.check_reference_consistency(genome_records, gff_records).ok


def test_feature_beyond_contig_end_is_an_error(genome_file, annotation_file):
    with annotation_file.open("a") as handle:
        handle.write("pPlasmid1\tRefSeq\tCDS\t1900\t2500\t.\t+\t0\tID=oops;locus_tag=b9999\n")

    _, genome_records = checks.check_genome(genome_file)
    _, gff_records = checks.check_annotation(annotation_file)
    report = checks.check_reference_consistency(genome_records, gff_records)

    assert "REFERENCE_COORDINATE_OUT_OF_BOUNDS" in check_ids(report, Severity.ERROR)


def test_genome_sequence_without_annotation_warns(genome_file, annotation_file):
    kept = [
        line
        for line in annotation_file.read_text().splitlines()
        if not line.startswith("pPlasmid1")
    ]
    annotation_file.write_text("\n".join(kept) + "\n")

    _, genome_records = checks.check_genome(genome_file)
    _, gff_records = checks.check_annotation(annotation_file)
    report = checks.check_reference_consistency(genome_records, gff_records)

    assert report.ok
    assert "REFERENCE_UNANNOTATED_SEQUENCE" in check_ids(report, Severity.WARNING)


# --------------------------------------------------------------------------
# Annotation structure
# --------------------------------------------------------------------------


def test_embedded_fasta_section_is_rejected(annotation_file):
    with annotation_file.open("a") as handle:
        handle.write("##FASTA\n>NC_000913.3\nACGTACGT\n")

    report, _ = checks.check_annotation(annotation_file)
    assert "ANNOTATION_EMBEDDED_FASTA" in check_ids(report, Severity.ERROR)


def test_embedded_fasta_does_not_break_parsing(annotation_file):
    """The FASTA block must not be parsed as if it were feature rows."""
    original, _, _ = parse_annotation(annotation_file)
    with annotation_file.open("a") as handle:
        handle.write("##FASTA\n>NC_000913.3\nACGTACGT\n")

    records, _, has_fasta = parse_annotation(annotation_file)
    assert has_fasta
    assert len(records) == len(original)


def test_annotation_without_cds_is_an_error(annotation_file):
    annotation_file.write_text(annotation_file.read_text().replace("\tCDS\t", "\texon\t"))
    report, _ = checks.check_annotation(annotation_file)
    assert "ANNOTATION_NO_CDS" in check_ids(report, Severity.ERROR)


def test_annotation_without_structural_rna_warns(annotation_file):
    text = annotation_file.read_text().replace("\trRNA\t", "\tCDS\t").replace("\ttRNA\t", "\tCDS\t")
    annotation_file.write_text(text)
    report, _ = checks.check_annotation(annotation_file)
    assert "ANNOTATION_NO_RRNA" in check_ids(report, Severity.WARNING)


def test_invalid_strand_is_an_error(annotation_file):
    annotation_file.write_text(annotation_file.read_text().replace("\t.\t+\t", "\t.\tX\t", 1))
    report, _ = checks.check_annotation(annotation_file)
    assert "ANNOTATION_BAD_STRAND" in check_ids(report, Severity.ERROR)


# --------------------------------------------------------------------------
# Genome file
# --------------------------------------------------------------------------


def test_duplicate_genome_identifiers_are_rejected(genome_file):
    genome_file.write_text(genome_file.read_text() + ">NC_000913.3 again\nACGTACGT\n")
    report, _ = checks.check_genome(genome_file)
    assert "GENOME_DUPLICATE_IDS" in check_ids(report, Severity.ERROR)


def test_protein_fasta_is_rejected(tmp_path: Path):
    path = tmp_path / "proteins.fa"
    path.write_text(">p1 some protein\nMKQLEDKVEELLSKNYHLENEVARLKKLV\n")
    report, _ = checks.check_genome(path)
    assert "GENOME_INVALID_CHARACTERS" in check_ids(report, Severity.ERROR)


def test_missing_genome_is_reported_not_raised(tmp_path: Path):
    report, records = checks.check_genome(tmp_path / "absent.fa")
    assert "GENOME_MISSING" in check_ids(report, Severity.ERROR)
    assert records == []


def test_fasta_identifier_is_first_token_only(genome_file):
    records = parse_fasta(genome_file)
    assert [r.identifier for r in records] == ["NC_000913.3", "pPlasmid1"]
    assert records[0].header.startswith("NC_000913.3 Escherichia")


# --------------------------------------------------------------------------
# Sample sheet and fastq files
# --------------------------------------------------------------------------


def test_duplicate_sample_key_is_an_error(samples):
    samples.loc[samples.index[1], "replicate"] = samples.loc[samples.index[0], "replicate"]
    report = checks.check_sample_sheet(samples)
    assert "SAMPLES_DUPLICATE_KEY" in check_ids(report, Severity.ERROR)


def test_hyphen_in_condition_is_an_error(samples):
    samples.loc[samples.index[0], "condition"] = "A-x"
    report = checks.check_sample_sheet(samples)
    assert "SAMPLES_CONDITION_NOT_ALNUM" in check_ids(report, Severity.ERROR)


def test_reused_fastq_file_warns(samples):
    samples.loc[samples.index[1], "fastqFile"] = samples.loc[samples.index[0], "fastqFile"]
    report = checks.check_sample_sheet(samples)
    assert "SAMPLES_REUSED_FASTQ" in check_ids(report, Severity.WARNING)


def test_missing_fastq_is_an_error(samples):
    samples.loc[samples.index[0], "fastqFile"] = "fastq/does-not-exist.fastq.gz"
    report = checks.check_fastq_files(samples)
    assert "FASTQ_MISSING" in check_ids(report, Severity.ERROR)


def test_truncated_gzip_is_an_error(samples):
    path = Path(samples.loc[samples.index[0], "fastqFile"])
    path.write_bytes(path.read_bytes()[:12])
    report = checks.check_fastq_files(samples)
    assert "FASTQ_UNREADABLE" in check_ids(report, Severity.ERROR)


# --------------------------------------------------------------------------
# Config semantics
# --------------------------------------------------------------------------


def test_differential_expression_without_rna_is_an_error(config, samples):
    config["differentialExpressionSettings"]["differentialExpression"] = "on"
    ribo_only = samples[samples["method"] == "RIBO"]
    report = checks.check_config_semantics(config, ribo_only)
    assert "DIFFEX_NO_RNA" in check_ids(report, Severity.ERROR)


def test_contrast_naming_unknown_condition_is_an_error(config, samples):
    config["differentialExpressionSettings"]["differentialExpression"] = "on"
    config["differentialExpressionSettings"]["contrasts"] = ["B-Z"]
    report = checks.check_config_semantics(config, samples)
    assert "DIFFEX_UNKNOWN_CONDITION" in check_ids(report, Severity.ERROR)


def test_single_condition_with_diffex_is_an_error(config, samples):
    config["differentialExpressionSettings"]["differentialExpression"] = "on"
    report = checks.check_config_semantics(config, samples[samples["condition"] == "A"])
    assert "DIFFEX_SINGLE_CONDITION" in check_ids(report, Severity.ERROR)


def test_single_replicate_warns(config, samples):
    config["differentialExpressionSettings"]["differentialExpression"] = "on"
    reduced = samples[samples["replicate"] == "1"]
    report = checks.check_config_semantics(config, reduced)
    assert "DIFFEX_SINGLE_REPLICATE" in check_ids(report, Severity.WARNING)


def test_invalid_adapter_is_an_error(config, samples):
    config["biologySettings"]["adapterS3"] = "AGATCGGXYZ"
    report = checks.check_config_semantics(config, samples)
    assert "CONFIG_INVALID_ADAPTER" in check_ids(report, Severity.ERROR)


def test_metagene_window_longer_than_every_gene_is_an_error(config, annotation_file):
    config["metageneSettings"]["positionsInORF"] = 5000
    _, gff_records = checks.check_annotation(annotation_file)
    report = checks.check_metagene_annotation_coverage(config, gff_records)
    assert "METAGENE_NO_GENES_SURVIVE" in check_ids(report, Severity.ERROR)


# --------------------------------------------------------------------------
# Reporting behaviour
# --------------------------------------------------------------------------


def test_all_errors_are_reported_at_once(config, samples):
    """A run must not stop at the first problem."""
    config["differentialExpressionSettings"]["differentialExpression"] = "on"
    config["differentialExpressionSettings"]["contrasts"] = ["B-Z"]
    config["biologySettings"]["adapterS3"] = "NOTDNA!"
    report = checks.check_config_semantics(config, samples[samples["method"] == "RIBO"])

    assert {"DIFFEX_NO_RNA", "DIFFEX_UNKNOWN_CONDITION", "CONFIG_INVALID_ADAPTER"} <= check_ids(
        report, Severity.ERROR
    )


def test_raise_on_error_mentions_every_error(config, samples):
    config["biologySettings"]["adapterS3"] = "NOTDNA!"
    report = checks.check_config_semantics(config, samples)
    with pytest.raises(ValueError, match="CONFIG_INVALID_ADAPTER"):
        report.raise_on_error()


def test_raise_on_error_is_silent_when_clean(config, samples):
    checks.check_config_semantics(config, samples).raise_on_error()
