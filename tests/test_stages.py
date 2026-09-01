"""Tests for the stage selection that decides which parts of the pipeline run."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

from lib import checks
from lib.stages import (
    PRESETS,
    RIBO_LIKE_STAGES,
    STAGE_NAMES,
    StageError,
    resolve_stages,
)
from lib.validation import Severity


def check_ids(report, severity=None):
    return {
        f.check for f in report.findings if severity is None or f.severity is severity
    }


def test_full_preset_selects_every_stage():
    assert resolve_stages({"workflowSettings": {"stages": "full"}}) == STAGE_NAMES


def test_preprocessing_preset_stops_at_the_qc_report():
    assert resolve_stages({"workflowSettings": {"stages": "preprocessing"}}) == [
        "trimming",
        "mapping",
        "qc",
    ]


def test_a_single_stage_is_not_expanded_into_its_preset():
    """"mapping" names a stage, so it must not drag in the trimming reports."""
    assert resolve_stages({"workflowSettings": {"stages": ["mapping"]}}) == ["mapping"]


def test_stages_are_returned_in_pipeline_order():
    requested = ["overview", "qc", "mapping"]
    assert resolve_stages({"workflowSettings": {"stages": requested}}) == [
        "mapping",
        "qc",
        "overview",
    ]


def test_command_line_override_wins_over_the_config_file():
    config = {"stages": "mapping,tracks", "workflowSettings": {"stages": "full"}}
    assert resolve_stages(config) == ["mapping", "tracks"]


def test_ribo_only_stages_are_dropped_without_ribo_libraries():
    config = {"workflowSettings": {"stages": ["mapping", "predictions", "overview"]}}
    assert resolve_stages(config, has_ribo=False) == ["mapping"]


@pytest.mark.parametrize(
    "selection",
    ["full", ["mapping", "metagene", "tis_advisor"]],
    ids=["default-full", "explicit"],
)
def test_rna_only_projects_drop_metagene_stages(selection):
    stages = resolve_stages(
        {"workflowSettings": {"stages": selection}},
        has_ribo=False,
        has_ribo_like=False,
    )

    assert not RIBO_LIKE_STAGES.intersection(stages)
    assert "mapping" in stages


def test_tis_libraries_keep_metagene_stages_without_general_ribo():
    config = {"workflowSettings": {"stages": ["metagene", "tis_advisor"]}}

    assert resolve_stages(
        config,
        has_ribo=False,
        has_ribo_like=True,
    ) == ["metagene", "tis_advisor"]


def test_missing_setting_defaults_to_a_full_run():
    assert resolve_stages({}) == STAGE_NAMES


@pytest.mark.parametrize("value", ["bam_files", ["mapping", "bam_files"], "", 12])
def test_unusable_stage_selections_are_rejected(value):
    with pytest.raises(StageError):
        resolve_stages({"workflowSettings": {"stages": value}})


def test_presets_only_contain_real_stages():
    for stages in PRESETS.values():
        assert set(stages) <= set(STAGE_NAMES)


def test_schema_and_module_agree_on_the_stage_names():
    """The schema rejects unknown stages, so its enum must not drift."""
    schema = yaml.safe_load(
        (Path(__file__).resolve().parent.parent / "workflow" / "schemas" / "config.schema.yaml").read_text()
    )
    assert schema["$defs"]["stageName"]["enum"] == STAGE_NAMES


def test_the_shipped_config_selects_usable_stages():
    config = yaml.safe_load(
        (Path(__file__).resolve().parent.parent / "config" / "config.yaml").read_text()
    )
    assert resolve_stages(config)


# --------------------------------------------------------------------------
# How the preflight reacts
# --------------------------------------------------------------------------


def test_unknown_stage_is_reported_by_the_preflight(config, samples):
    config["workflowSettings"]["stages"] = ["bam_files"]
    report = checks.check_config_semantics(config, samples)
    assert "CONFIG_UNKNOWN_STAGE" in check_ids(report, Severity.ERROR)


def test_predictions_without_ribo_libraries_warn(config, samples):
    config["workflowSettings"]["stages"] = ["mapping", "predictions"]
    rna_only = samples[samples["method"] != "RIBO"]
    report = checks.check_config_semantics(config, rna_only)
    assert "STAGES_NEED_RIBO" in check_ids(report, Severity.WARNING)


def test_rna_only_metagene_stages_warn_and_are_skipped(config, samples):
    config["workflowSettings"]["stages"] = ["mapping", "metagene", "tis_advisor"]
    rna_only = samples[samples["method"] == "RNA"]
    report = checks.check_config_semantics(config, rna_only)

    findings = [
        finding
        for finding in report.findings
        if finding.check == "STAGES_NEED_RIBO"
    ]
    assert len(findings) == 1
    assert set(findings[0].items) == {"metagene", "tis_advisor"}


def test_rna_only_skips_metagene_annotation_coverage(config, annotation_file):
    config["workflowSettings"]["stages"] = "full"
    config["metageneSettings"]["positionsInORF"] = 5000
    _, records = checks.check_annotation(annotation_file)

    report = checks.check_metagene_annotation_coverage(
        config,
        records,
        methods={"RNA"},
    )

    assert "METAGENE_NO_GENES_SURVIVE" not in check_ids(report, Severity.ERROR)


def test_diffex_settings_are_not_checked_when_the_stage_is_off(config, samples):
    """A contrast naming an unknown condition is harmless if diffex never runs."""
    config["workflowSettings"]["stages"] = ["mapping"]
    config["differentialExpressionSettings"]["contrasts"] = ["B-Z"]
    report = checks.check_config_semantics(config, samples)
    assert "DIFFEX_UNKNOWN_CONDITION" not in check_ids(report)


@pytest.fixture(scope="module")
def snakemake_command():
    executable = shutil.which("snakemake")
    if executable:
        return [executable]

    conda = shutil.which("conda")
    if conda:
        command = [conda, "run", "-n", "snakemake", "snakemake"]
        probe = subprocess.run(
            [*command, "--version"],
            capture_output=True,
            text=True,
        )
        if probe.returncode == 0:
            return command

    pytest.skip("Snakemake is not available for the RNA-only DAG test")


@pytest.mark.parametrize(
    "selection",
    ["full", ["mapping", "metagene", "tis_advisor"]],
    ids=["default-full", "explicit"],
)
def test_rna_only_dag_omits_metagene_commands(
    selection,
    snakemake_command,
    genome_file,
    annotation_file,
    samples,
    tmp_path,
):
    """An RNA-only dry-run must never render an empty metagene ``-a`` command."""
    repo = Path(__file__).resolve().parent.parent
    workflow_config = yaml.safe_load((repo / "config" / "config.yaml").read_text())

    sample_path = tmp_path / "samples.tsv"
    samples[samples["method"] == "RNA"].fillna("").to_csv(
        sample_path,
        sep="\t",
        index=False,
    )
    workflow_config["biologySettings"].update(
        {
            "genome": str(genome_file),
            "annotation": str(annotation_file),
            "samples": str(sample_path),
        }
    )
    workflow_config["workflowSettings"]["stages"] = selection

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(workflow_config, sort_keys=False))

    result = subprocess.run(
        [
            *snakemake_command,
            "--dry-run",
            "--printshellcmds",
            "--cores",
            "1",
            "--snakefile",
            str(repo / "workflow" / "Snakefile"),
            "--directory",
            str(tmp_path),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / ".cache")},
    )
    rendered_dag = result.stdout + result.stderr

    assert result.returncode == 0, rendered_dag
    assert "readLengthStatistics" not in rendered_dag
    assert "read_length_statistics.py" not in rendered_dag
    assert "metageneProfiling" not in rendered_dag
    assert "metagene_profiling.py" not in rendered_dag
    assert "tisAdvisor" not in rendered_dag
    assert "tis_advisor.py" not in rendered_dag
