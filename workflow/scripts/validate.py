"""
Entry points for validating the HRIBO configuration and inputs.

Structural validation of the config and the sample sheet is delegated to the
JSON schemas in workflow/schemas; everything that a schema cannot express --
consistency between the genome and the annotation, whether the requested
contrasts are achievable, whether the fastq files are readable -- lives in
lib.checks.

All of this runs at parse time, before any job is scheduled, so that bad input
fails immediately and with an explanation rather than several hours later inside
an unrelated tool.

Author: Rick Gelhausen
"""

from __future__ import annotations

import sys
from pathlib import Path

from lib import checks
from lib.stages import (
    ANNOTATION_REQUIRED_STAGES,
    FASTQ_REQUIRED_STAGES,
    GENOME_REQUIRED_STAGES,
    RIBO_LIKE_METHODS,
    StageError,
    resolve_stages,
)
from lib.validation import Severity, ValidationReport

SCHEMA_DIR = Path(__file__).resolve().parent.parent / "schemas"


def validate_config(config, unique_conditions=None):
    """Validate the config file against the JSON schema.

    `unique_conditions` is accepted for backwards compatibility with the
    previous signature; contrast-versus-condition checking now happens in
    validate_inputs, where the sample sheet is available.
    """
    from snakemake.utils import validate as schema_validate

    schema_validate(config, str(SCHEMA_DIR / "config.schema.yaml"))
    print("Config file validated against schema.", file=sys.stderr)


def validate_sample_sheet(samples):
    """Validate the sample sheet against the JSON schema, plus cross-row checks."""
    from snakemake.utils import validate as schema_validate

    # Snakemake removes pandas NULL values before validating each row.
    # Consequently, a nullable field alone cannot distinguish an absent header
    # from an intentionally empty single-end cell. Check both path headers
    # before normalization so a malformed sheet cannot be concealed.
    for column in ("fastqFile", "fastqFile2"):
        if column not in samples.columns:
            guidance = (
                " Keep the column and leave its cells empty for single-end "
                "libraries."
                if column == "fastqFile2"
                else ""
            )
            raise ValueError(
                f"Sample sheet validation: missing required column {column!r}."
                f"{guidance}"
            )

    # Preserve the required headers while representing pandas' several NULL
    # spellings as the empty string understood by all downstream layout checks.
    # Strip only surrounding whitespace: spaces inside a path remain valid, but
    # preflight and DAG construction must see exactly the same pathname.
    for column in ("fastqFile", "fastqFile2"):
        samples[column] = samples[column].astype("object").where(
            samples[column].notna(), ""
        )
        samples[column] = samples[column].map(
            lambda value: value.strip() if isinstance(value, str) else value
        )

    schema_validate(samples, str(SCHEMA_DIR / "samples.schema.yaml"))

    report = checks.check_sample_sheet(samples)
    report.raise_on_error("Sample sheet validation")
    _emit(report)
    print("Sample sheet validated.", file=sys.stderr)


def validate_inputs(config, samples, verbose: bool = True) -> ValidationReport:
    """Validate only the inputs required by the requested workflow stages.

    Returns the report so that callers can render it; raises on any error.
    """
    biology = config.get("biologySettings", {})
    genome_path = Path(biology.get("genome", ""))
    annotation_path = Path(biology.get("annotation", ""))

    report = ValidationReport()

    methods = set(samples["method"].astype(str))
    try:
        stages = resolve_stages(
            config,
            has_ribo="RIBO" in methods,
            has_ribo_like=bool(RIBO_LIKE_METHODS & methods),
        )
    except StageError:
        # check_config_semantics records the actionable stage error below.
        stages = []

    needs_genome = bool(GENOME_REQUIRED_STAGES.intersection(stages))
    needs_annotation = bool(ANNOTATION_REQUIRED_STAGES.intersection(stages))
    needs_fastq = bool(FASTQ_REQUIRED_STAGES.intersection(stages))

    genome_records = []
    if needs_genome:
        genome_report, genome_records = checks.check_genome(genome_path)
        report.extend(genome_report)

    prediction_settings = config.get("predictionSettings", {})
    deepribo_setting = prediction_settings.get("deepribo", "")
    runs_deepribo = (
        isinstance(deepribo_setting, str)
        and deepribo_setting.lower() == "on"
        and bool({"predictions", "overview"}.intersection(stages))
    )
    if runs_deepribo and genome_records:
        report.extend(checks.check_deepribo_genome(genome_records))

    gff_records = []
    if needs_annotation:
        annotation_report, gff_records = checks.check_annotation(annotation_path)
        report.extend(annotation_report)

    if needs_genome and needs_annotation:
        report.extend(checks.check_reference_consistency(genome_records, gff_records))
    if needs_fastq:
        report.extend(checks.check_fastq_files(samples))
    report.extend(checks.check_config_semantics(config, samples))
    if needs_annotation:
        report.extend(
            checks.check_metagene_annotation_coverage(
                config,
                gff_records,
                methods=methods,
            )
        )

    if verbose:
        _emit(report)

    report.raise_on_error("Input validation")
    return report


def _emit(report: ValidationReport) -> None:
    """Print warnings and informational findings; errors are raised separately."""
    for finding in report.sorted_findings():
        if finding.severity is not Severity.ERROR:
            # stderr, so that machine-readable stdout (--summary, --list) stays clean
            print(finding.to_text(), file=sys.stderr)
