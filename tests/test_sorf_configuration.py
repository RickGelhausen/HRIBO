"""Configuration compatibility and preflight checks for separate sORF profiles."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import jsonschema
import pytest
import yaml

from lib import checks
from lib.validation import GffRecord, Severity


REPO = Path(__file__).resolve().parent.parent


def _short_cds():
    return GffRecord(
        line_number=2,
        seqid="chr1",
        source="test",
        feature="CDS",
        start=101,
        end=160,
        score=".",
        strand="+",
        phase="0",
        attributes="ID=short",
    )


def _length_finding(report):
    return next(
        finding
        for finding in report.findings
        if finding.check == "METAGENE_NO_GENES_SURVIVE"
    )


@pytest.mark.parametrize("limit", [None, 0, 300, 450])
def test_sorf_schema_accepts_old_and_enabled_configurations(limit):
    workflow_config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    schema = yaml.safe_load(
        (REPO / "workflow" / "schemas" / "config.schema.yaml").read_text()
    )
    if limit is None:
        workflow_config["metageneSettings"].pop("sorfMaxLength")
    else:
        workflow_config["metageneSettings"]["sorfMaxLength"] = limit

    jsonschema.Draft202012Validator(schema).validate(workflow_config)


@pytest.mark.parametrize("limit", [-1, 299.5, "300", True])
def test_sorf_schema_rejects_invalid_limits(limit):
    workflow_config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    schema = yaml.safe_load(
        (REPO / "workflow" / "schemas" / "config.schema.yaml").read_text()
    )
    workflow_config["metageneSettings"]["sorfMaxLength"] = limit

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(workflow_config)


@pytest.mark.parametrize("stages", [["metagene"], ["metagene", "tis_advisor"]])
def test_sorf_group_allows_an_annotation_of_only_short_cdss(config, stages):
    config["workflowSettings"]["stages"] = stages
    config["metageneSettings"]["sorfMaxLength"] = 300

    report = checks.check_metagene_annotation_coverage(config, [_short_cds()])

    assert report.ok
    finding = _length_finding(report)
    assert finding.severity is Severity.WARNING
    if "tis_advisor" in stages:
        assert "advisor" in finding.detail


@pytest.mark.parametrize("limit", [None, 0, 60])
def test_no_usable_sorf_group_keeps_the_original_empty_annotation_error(
    config, limit
):
    config["workflowSettings"]["stages"] = ["metagene"]
    if limit is None:
        config["metageneSettings"].pop("sorfMaxLength", None)
    else:
        config["metageneSettings"]["sorfMaxLength"] = limit

    report = checks.check_metagene_annotation_coverage(config, [_short_cds()])

    assert _length_finding(report).severity is Severity.ERROR


def test_sorf_setting_does_not_relax_advisor_only_length_checks(config):
    config["workflowSettings"]["stages"] = ["tis_advisor"]
    config["metageneSettings"]["sorfMaxLength"] = 300

    report = checks.check_metagene_annotation_coverage(config, [_short_cds()])

    assert _length_finding(report).severity is Severity.ERROR


def test_length_preflight_is_inactive_when_the_length_filter_is_disabled(
    config, samples
):
    config["workflowSettings"]["stages"] = ["metagene", "tis_advisor"]
    config["metageneSettings"]["filteringMethods"] = ["rpkm"]

    coverage = checks.check_metagene_annotation_coverage(config, [_short_cds()])
    semantics = checks.check_config_semantics(config, samples)

    assert not coverage.findings
    assert "METAGENE_LENGTH_CUTOFF" not in {
        finding.check for finding in semantics.findings
    }


def test_length_preflight_uses_the_configured_effective_minimum(config):
    config["workflowSettings"]["stages"] = ["metagene"]
    config["metageneSettings"].update({"positionsInORF": 30, "lengthCutoff": 90})

    report = checks.check_metagene_annotation_coverage(config, [_short_cds()])

    assert _length_finding(report).severity is Severity.ERROR


@pytest.mark.parametrize("limit", [None, 0, 450])
def test_workflow_passes_the_sorf_limit_only_to_metagene_profiling(
    limit, snakemake_command, genome_file, annotation_file, samples, tmp_path
):
    workflow_config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    sample_path = tmp_path / "samples.tsv"
    library = samples[samples["method"] == "RIBO"].iloc[[0]].copy()
    library.fillna("").to_csv(sample_path, sep="\t", index=False)
    workflow_config["biologySettings"].update(
        {
            "genome": str(genome_file),
            "annotation": str(annotation_file),
            "samples": str(sample_path),
        }
    )
    workflow_config["workflowSettings"]["stages"] = ["metagene", "tis_advisor"]
    if limit is None:
        workflow_config["metageneSettings"].pop("sorfMaxLength")
    else:
        workflow_config["metageneSettings"]["sorfMaxLength"] = limit
    config_path = tmp_path / "config" / "config.yaml"
    config_path.parent.mkdir()
    config_path.write_text(yaml.safe_dump(workflow_config, sort_keys=False))

    result = subprocess.run(
        [
            *snakemake_command,
            "--dry-run",
            "--printshellcmds",
            "--cores",
            "1",
            "--snakefile",
            str(REPO / "workflow" / "Snakefile"),
            "--directory",
            str(tmp_path),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / ".cache")},
        timeout=120,
    )
    rendered = result.stdout + result.stderr

    assert result.returncode == 0, rendered
    expected_limit = 0 if limit is None else limit
    assert rendered.count(f"--sorf_max_length {expected_limit}") == 1
