"""Tests for the input preflight validation.

Each test asserts on the check identifier rather than on message text, so that
wording can be improved without breaking the suite.
"""

import gzip
import random
from pathlib import Path

import pandas as pd
import pytest
from snakemake_interface_common.exceptions import WorkflowError

from lib import checks
from lib.validation import Severity, parse_annotation, parse_fasta
from validate import validate_inputs, validate_sample_sheet


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


def test_trimming_only_ignores_unrelated_reference_files(config, samples, tmp_path):
    config["workflowSettings"]["stages"] = ["trimming"]
    config["biologySettings"]["genome"] = str(tmp_path / "missing-genome.fa")
    config["biologySettings"]["annotation"] = str(
        tmp_path / "missing-annotation.gff"
    )

    report = validate_inputs(config, samples, verbose=False)

    assert not {
        "GENOME_MISSING",
        "ANNOTATION_MISSING",
    }.intersection(check_ids(report, Severity.ERROR))


def test_trimming_only_still_requires_fastq_files(config, samples, tmp_path):
    config["workflowSettings"]["stages"] = ["trimming"]
    missing_reads = samples.copy()
    missing_reads["fastqFile"] = str(tmp_path / "missing.fastq.gz")

    with pytest.raises(ValueError, match="FASTQ_MISSING"):
        validate_inputs(config, missing_reads, verbose=False)


def test_genome_tracks_ignore_annotation_and_fastq_files(
    genome_file, config, samples, tmp_path
):
    config["workflowSettings"]["stages"] = ["genome_tracks"]
    config["biologySettings"]["genome"] = str(genome_file)
    config["biologySettings"]["annotation"] = str(
        tmp_path / "missing-annotation.gff"
    )
    missing_reads = samples.copy()
    missing_reads["fastqFile"] = str(tmp_path / "missing.fastq.gz")

    report = validate_inputs(config, missing_reads, verbose=False)

    assert not {
        "ANNOTATION_MISSING",
        "FASTQ_MISSING",
    }.intersection(check_ids(report, Severity.ERROR))


def test_genome_tracks_still_require_the_genome(config, samples, tmp_path):
    config["workflowSettings"]["stages"] = ["genome_tracks"]
    config["biologySettings"]["genome"] = str(tmp_path / "missing-genome.fa")

    with pytest.raises(ValueError, match="GENOME_MISSING"):
        validate_inputs(config, samples, verbose=False)


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
    original, _, _, original_malformed = parse_annotation(annotation_file)
    with annotation_file.open("a") as handle:
        handle.write("##FASTA\n>NC_000913.3\nACGTACGT\n")

    records, _, has_fasta, malformed = parse_annotation(annotation_file)
    assert has_fasta
    assert malformed == original_malformed == []
    assert len(records) == len(original)


def test_malformed_annotation_rows_are_blocking_and_reported(annotation_file):
    with annotation_file.open("a") as handle:
        handle.write("chr1\ttest\tCDS\t1\t9\t.\t+\t0\n")
        handle.write("chr1\ttest\tCDS\tstart\t9\t.\t+\t0\tID=bad;\n")

    report, records = checks.check_annotation(annotation_file)
    findings = [
        finding
        for finding in report.findings
        if finding.check == "ANNOTATION_MALFORMED_ROWS"
    ]

    assert len(findings) == 1
    assert findings[0].severity is Severity.ERROR
    assert any("expected 9" in item for item in findings[0].items)
    assert any("non-integer coordinates" in item for item in findings[0].items)
    assert records  # valid rows remain available for the other consistency checks


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


def _append_genome_symbol(genome_file: Path, symbol: str) -> None:
    with genome_file.open("a") as handle:
        handle.write(f"{symbol}\n")


@pytest.mark.parametrize("stage", ["predictions", "overview"])
def test_deepribo_rejects_broader_iupac_for_stages_that_run_it(
    genome_file, config, samples, stage
):
    _append_genome_symbol(genome_file, "R")
    config["workflowSettings"]["stages"] = [stage]

    with pytest.raises(ValueError) as error:
        validate_inputs(config, samples, verbose=False)

    message = str(error.value)
    assert "GENOME_DEEPRIBO_ALPHABET" in message
    assert "A/C/G/T/N" in message
    assert "pPlasmid1 (R)" in message
    assert "predictionSettings.deepribo to 'off'" in message


@pytest.mark.parametrize(
    ("stage", "deepribo"),
    [("mapping", "on"), ("predictions", "off")],
)
def test_broader_iupac_remains_valid_when_deepribo_will_not_run(
    genome_file, config, samples, stage, deepribo
):
    _append_genome_symbol(genome_file, "R")
    config["workflowSettings"]["stages"] = [stage]
    config["predictionSettings"]["deepribo"] = deepribo

    report = validate_inputs(config, samples, verbose=False)

    assert "GENOME_DEEPRIBO_ALPHABET" not in check_ids(report, Severity.ERROR)
    assert "GENOME_INVALID_CHARACTERS" not in check_ids(report, Severity.ERROR)


def test_deepribo_rejects_lowercase_genome_sequence(genome_file, config, samples):
    lines = genome_file.read_text().splitlines()
    genome_file.write_text(
        "\n".join(line if line.startswith(">") else line.lower() for line in lines)
        + "\nn\n"
    )
    config["workflowSettings"]["stages"] = ["predictions"]

    with pytest.raises(ValueError) as error:
        validate_inputs(config, samples, verbose=False)

    message = str(error.value)
    assert "GENOME_DEEPRIBO_ALPHABET" in message
    assert "uppercase A/C/G/T/N" in message
    assert "Convert lowercase sequence to uppercase" in message


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


@pytest.mark.parametrize("column", ["fastqFile", "fastqFile2"])
def test_sample_sheet_requires_fastq_path_headers(samples, column):
    malformed = samples.drop(columns=column)

    with pytest.raises(ValueError, match=rf"missing required column '{column}'"):
        validate_sample_sheet(malformed)


@pytest.mark.parametrize("empty_value", [None, pd.NA, float("nan"), ""])
def test_sample_sheet_accepts_empty_fastq_file_2_cells(samples, empty_value):
    single_end = samples.copy()
    single_end["fastqFile2"] = empty_value

    validate_sample_sheet(single_end)

    # Validation gives every downstream DAG consumer one stable empty value.
    assert single_end["fastqFile2"].tolist() == [""] * len(single_end)


def test_sample_sheet_rejects_a_whitespace_only_primary_fastq(samples):
    malformed = samples.iloc[[0]].copy()
    malformed.loc[:, "fastqFile"] = " \t "

    with pytest.raises(WorkflowError, match="Error validating row 0"):
        validate_sample_sheet(malformed)


def test_sample_sheet_normalizes_outer_fastq_path_whitespace(samples):
    padded = samples.iloc[[0]].copy()
    original_first = padded.iloc[0]["fastqFile"]
    padded.loc[:, "fastqFile"] = f"  {original_first}\t"
    padded.loc[:, "fastqFile2"] = "  reads with spaces/mate 2.fastq.gz  "

    validate_sample_sheet(padded)

    assert padded.iloc[0]["fastqFile"] == original_first
    assert padded.iloc[0]["fastqFile2"] == "reads with spaces/mate 2.fastq.gz"


@pytest.mark.parametrize("replicate", ["1", "9", "10", "123456"])
def test_sample_sheet_accepts_canonical_positive_replicates(samples, replicate):
    one_sample = samples.iloc[[0]].copy()
    one_sample.loc[:, "replicate"] = replicate

    validate_sample_sheet(one_sample)


@pytest.mark.parametrize(
    "replicate",
    ["0", "00", "01", "+1", "-1", "1.0", " 1", "1 "],
)
def test_sample_sheet_rejects_noncanonical_replicates(samples, replicate):
    one_sample = samples.iloc[[0]].copy()
    one_sample.loc[:, "replicate"] = replicate

    with pytest.raises(WorkflowError, match="Error validating row 0"):
        validate_sample_sheet(one_sample)


def test_truncated_gzip_is_an_error(samples):
    path = Path(samples.loc[samples.index[0], "fastqFile"])
    path.write_bytes(path.read_bytes()[:12])
    report = checks.check_fastq_files(samples)
    assert "FASTQ_UNREADABLE" in check_ids(report, Severity.ERROR)


def _large_fastq_payload(record_count: int = 1024) -> bytes:
    """Return enough poorly-compressible FASTQ data to exceed gzip's read buffer."""
    rng = random.Random(0)
    records = []
    for index in range(record_count):
        sequence = "".join(rng.choice("ACGT") for _ in range(100))
        records.append(f"@read-{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    return "".join(records).encode("ascii")


def _samples_for_fastq(samples, path: Path):
    selected = samples.iloc[[0]].copy()
    selected.loc[:, "fastqFile"] = str(path)
    selected.loc[:, "fastqFile2"] = ""
    return selected


@pytest.mark.parametrize("suffix", [".fastq", ".fastq.gz"])
def test_gzip_content_is_accepted_regardless_of_source_suffix(samples, tmp_path, suffix):
    path = tmp_path / f"complete{suffix}"
    payload = _large_fastq_payload()
    path.write_bytes(gzip.compress(payload, mtime=0))

    report = checks.check_fastq_files(_samples_for_fastq(samples, path))

    assert report.ok, report.to_text()


@pytest.mark.parametrize("suffix", [".fastq", ".fastq.gz"])
def test_plain_fastq_content_is_rejected_regardless_of_suffix(
    samples, tmp_path, suffix
):
    path = tmp_path / f"plain{suffix}"
    path.write_bytes(b"@read\nACGT\n+\nIIII\n")

    report = checks.check_fastq_files(_samples_for_fastq(samples, path))

    assert "FASTQ_NOT_GZIP" in check_ids(report, Severity.ERROR)
    finding = next(f for f in report.findings if f.check == "FASTQ_NOT_GZIP")
    assert str(path) in finding.items


def test_reused_fastq_path_is_streamed_only_once(samples, monkeypatch):
    path = Path(samples.loc[samples.index[0], "fastqFile"])
    reused = samples.copy()
    reused.loc[:, "fastqFile"] = str(path)
    reused.loc[:, "fastqFile2"] = ""
    scanned = []
    validate_fastq = checks._fastq_format_problem

    def record_scan(candidate):
        scanned.append(candidate)
        return validate_fastq(candidate)

    monkeypatch.setattr(checks, "_fastq_format_problem", record_scan)

    report = checks.check_fastq_files(reused)

    assert report.ok, report.to_text()
    assert scanned == [path]


def test_gzip_truncated_after_the_initial_probe_is_an_error(samples, tmp_path):
    path = tmp_path / "late-truncation.fastq.gz"
    payload = _large_fastq_payload()
    compressed = gzip.compress(payload, mtime=0)
    assert len(payload) > 1024
    assert len(compressed) > 8192
    path.write_bytes(compressed[:-8])

    report = checks.check_fastq_files(_samples_for_fastq(samples, path))

    assert "FASTQ_UNREADABLE" in check_ids(report, Severity.ERROR)
    finding = next(f for f in report.findings if f.check == "FASTQ_UNREADABLE")
    assert str(path) in finding.items[0]


def test_gzip_crc_corruption_after_the_initial_probe_is_an_error(samples, tmp_path):
    path = tmp_path / "late-crc-corruption.fastq.gz"
    payload = _large_fastq_payload()
    compressed = bytearray(gzip.compress(payload, mtime=0))
    assert len(compressed) > 8192
    compressed[-8] ^= 1  # The first byte of the gzip trailer's CRC32.
    path.write_bytes(compressed)

    report = checks.check_fastq_files(_samples_for_fastq(samples, path))

    assert "FASTQ_UNREADABLE" in check_ids(report, Severity.ERROR)
    finding = next(f for f in report.findings if f.check == "FASTQ_UNREADABLE")
    assert "CRC check failed" in finding.items[0]


@pytest.mark.parametrize("suffix", [".fastq", ".fastq.gz"])
@pytest.mark.parametrize(
    ("malformed_record", "reason"),
    [
        pytest.param(b"read-81\nACGT\n+\nIIII\n", "header", id="header-marker"),
        pytest.param(b"@read-81\nACGT\nseparator\nIIII\n", "separator", id="plus-marker"),
        pytest.param(b"@read-81\nACGT\n+\nIII\n", "quality", id="length-mismatch"),
        pytest.param(b"@read-81\nACGT\n+\n", "incomplete", id="incomplete-record"),
    ],
)
def test_malformed_record_after_the_initial_probe_is_an_error(
    samples, tmp_path, suffix, malformed_record, reason
):
    path = tmp_path / f"malformed-late{suffix}"
    prefix = b"@valid\nACGTACGT\n+\nIIIIIIII\n" * 80
    assert len(prefix) > 1024
    payload = prefix + malformed_record
    path.write_bytes(gzip.compress(payload, mtime=0))

    report = checks.check_fastq_files(_samples_for_fastq(samples, path))

    assert "FASTQ_MALFORMED" in check_ids(report, Severity.ERROR)
    finding = next(f for f in report.findings if f.check == "FASTQ_MALFORMED")
    assert str(path) in finding.items[0]
    assert "record 81" in finding.items[0]
    assert reason in finding.items[0]


# --------------------------------------------------------------------------
# Config semantics
# --------------------------------------------------------------------------


def test_differential_expression_without_rna_is_an_error(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    ribo_only = samples[samples["method"] == "RIBO"]
    report = checks.check_config_semantics(config, ribo_only)
    assert "DIFFEX_NO_RNA" in check_ids(report, Severity.ERROR)


def test_contrast_naming_unknown_condition_is_an_error(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    config["differentialExpressionSettings"]["contrasts"] = ["B-Z"]
    report = checks.check_config_semantics(config, samples)
    assert "DIFFEX_UNKNOWN_CONDITION" in check_ids(report, Severity.ERROR)


def test_self_contrast_is_an_error(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    config["differentialExpressionSettings"]["contrasts"] = ["A-A"]

    report = checks.check_config_semantics(config, samples)

    assert "DIFFEX_SELF_CONTRAST" in check_ids(report, Severity.ERROR)


def test_single_condition_with_diffex_is_an_error(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    report = checks.check_config_semantics(config, samples[samples["condition"] == "A"])
    assert "DIFFEX_SINGLE_CONDITION" in check_ids(report, Severity.ERROR)


def test_single_replicate_diffex_is_an_error(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    reduced = samples[samples["replicate"] == "1"]
    report = checks.check_config_semantics(config, reduced)
    assert "DIFFEX_SINGLE_REPLICATE" in check_ids(report, Severity.ERROR)


def test_unselected_single_replicate_condition_does_not_block_diffex(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    config["differentialExpressionSettings"]["contrasts"] = ["B-A"]
    extra = samples[
        (samples["condition"] == "A") & (samples["replicate"] == "1")
    ].copy()
    extra["condition"] = "C"

    report = checks.check_config_semantics(
        config, pd.concat([samples, extra], ignore_index=True)
    )

    assert "DIFFEX_SINGLE_REPLICATE" not in check_ids(report, Severity.ERROR)


def test_auto_contrasts_require_replicates_for_every_condition(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    extra = samples[
        (samples["condition"] == "A") & (samples["replicate"] == "1")
    ].copy()
    extra["condition"] = "C"

    report = checks.check_config_semantics(
        config, pd.concat([samples, extra], ignore_index=True)
    )

    assert "DIFFEX_SINGLE_REPLICATE" in check_ids(report, Severity.ERROR)


def test_auto_contrasts_ignore_conditions_without_matched_ribo_and_rna(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    tis_only = samples[
        (samples["method"] == "RIBO") & (samples["condition"] == "A")
    ].copy()
    tis_only["method"] = "TIS"
    tis_only["condition"] = "TIS-only"

    report = checks.check_config_semantics(
        config, pd.concat([samples, tis_only], ignore_index=True)
    )

    assert "DIFFEX_MISSING_GROUP" not in check_ids(report, Severity.ERROR)
    assert "DIFFEX_SINGLE_REPLICATE" not in check_ids(report, Severity.ERROR)


def test_diffex_rejects_unequal_ribo_and_rna_table_widths(config, samples):
    config["workflowSettings"]["stages"] = ["differential_expression"]
    config["differentialExpressionSettings"]["contrasts"] = ["B-A"]
    unbalanced = samples[
        ~(
            (samples["method"] == "RNA")
            & (samples["condition"] == "B")
            & (samples["replicate"] == "2")
        )
    ]

    report = checks.check_config_semantics(config, unbalanced)

    assert "DIFFEX_UNMATCHED_REPLICATES" in check_ids(report, Severity.ERROR)


def test_pca_stage_rejects_a_single_library(config, samples):
    config["workflowSettings"]["stages"] = ["pca"]
    report = checks.check_config_semantics(config, samples.iloc[[0]])

    assert "PCA_TOO_FEW_LIBRARIES" in check_ids(report, Severity.ERROR)


def test_pca_stage_accepts_two_libraries(config, samples):
    config["workflowSettings"]["stages"] = ["pca"]
    report = checks.check_config_semantics(config, samples.iloc[:2])

    assert "PCA_TOO_FEW_LIBRARIES" not in check_ids(report, Severity.ERROR)


def test_single_library_is_allowed_when_pca_is_not_selected(config, samples):
    config["workflowSettings"]["stages"] = ["mapping"]
    report = checks.check_config_semantics(config, samples.iloc[[0]])

    assert "PCA_TOO_FEW_LIBRARIES" not in check_ids(report, Severity.ERROR)


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
    config["workflowSettings"]["stages"] = ["differential_expression"]
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
