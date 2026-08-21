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

    schema_validate(samples, str(SCHEMA_DIR / "samples.schema.yaml"))

    report = checks.check_sample_sheet(samples)
    report.raise_on_error("Sample sheet validation")
    _emit(report)
    print("Sample sheet validated.", file=sys.stderr)


def validate_inputs(config, samples, verbose: bool = True) -> ValidationReport:
    """Run the full preflight over the reference files, fastq files and config.

    Returns the report so that callers can render it; raises on any error.
    """
    biology = config.get("biologySettings", {})
    genome_path = Path(biology.get("genome", ""))
    annotation_path = Path(biology.get("annotation", ""))

    report = ValidationReport()

    genome_report, genome_records = checks.check_genome(genome_path)
    report.extend(genome_report)

    annotation_report, gff_records = checks.check_annotation(annotation_path)
    report.extend(annotation_report)

    report.extend(checks.check_reference_consistency(genome_records, gff_records))
    report.extend(checks.check_fastq_files(samples))
    report.extend(checks.check_config_semantics(config, samples))
    report.extend(checks.check_metagene_annotation_coverage(config, gff_records))

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
