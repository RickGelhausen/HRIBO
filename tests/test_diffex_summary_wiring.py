"""The condition overview is a differential-expression target, not a prediction target."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import yaml


REPO = Path(__file__).resolve().parent.parent


def test_detection_thresholds_are_shipped_and_validated():
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    schema = yaml.safe_load(
        (REPO / "workflow" / "schemas" / "config.schema.yaml").read_text()
    )
    settings = config["differentialExpressionSettings"]
    settings_schema = schema["properties"]["differentialExpressionSettings"]

    for key, expected, kind, lower_bound in (
        ("detectionMinCPM", 1, "number", 0),
        ("detectionMinCount", 10, "integer", 1),
        ("detectionMinReplicates", 2, "integer", 1),
    ):
        field = settings_schema["properties"][key]
        assert settings[key] == expected
        assert field["type"] == kind
        assert field["minimum"] == lower_bound
        assert field["default"] == expected
        assert key not in settings_schema["required"]


def test_diffex_stage_requests_condition_overview_without_predictions(
    snakemake_command,
    genome_file,
    annotation_file,
    samples,
    tmp_path,
):
    workdir = tmp_path / "diffex stage with spaces"
    workdir.mkdir()
    sample_path = workdir / "samples.tsv"
    samples.fillna("").to_csv(sample_path, sep="\t", index=False)

    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["biologySettings"].update(
        {
            "genome": str(genome_file),
            "annotation": str(annotation_file),
            "samples": str(sample_path),
        }
    )
    config["workflowSettings"]["stages"] = ["differential_expression"]
    config["differentialExpressionSettings"]["contrasts"] = ["A-B"]
    config_path = workdir / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    result = subprocess.run(
        [
            *snakemake_command,
            "all",
            "--dry-run",
            "--printshellcmds",
            "--cores",
            "1",
            "--snakefile",
            str(REPO / "workflow" / "Snakefile"),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "XDG_CACHE_HOME": str(workdir / ".cache")},
        timeout=120,
    )
    rendered = result.stdout + result.stderr
    assert result.returncode == 0, rendered
    assert "rule conditionOverview:" in rendered
    assert "diffex_summary/condition_overview.html" in rendered
    assert "diffex_summary/condition_matrix.tsv" in rendered
    assert "diffex_summary/contrast_matrix.tsv" in rendered
    assert "diffex_summary/browser_tracks.tsv" in rendered
    assert "diffex_summary/browser" in rendered
    assert "--contrasts A-B" in rendered
    assert "--min_cpm 1" in rendered
    assert "--min_count 10" in rendered
    assert "--min_replicates 2" in rendered
    assert "rule poolxtail:" in rendered
    assert "rule poolriborex:" in rendered
    assert "rule pooldeltate:" in rendered
    assert "rule createOverviewTable:" not in rendered
    assert "rule asiteOccupancy:" not in rendered
    assert "rule reparation:" not in rendered


def test_condition_overview_rule_executes_from_existing_diffex_inputs(
    snakemake_command, genome_file, annotation_file, samples, tmp_path
):
    workdir = tmp_path / "summary rule with spaces"
    workdir.mkdir()
    sample_path = workdir / "samples.tsv"
    samples.fillna("").to_csv(sample_path, sep="\t", index=False)
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["biologySettings"].update({
        "genome": str(genome_file), "annotation": str(annotation_file),
        "samples": str(sample_path),
    })
    config["workflowSettings"]["stages"] = ["differential_expression"]
    config["differentialExpressionSettings"]["contrasts"] = ["A-B"]
    config_path = workdir / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    readcounts = workdir / "readcounts"
    readcounts.mkdir()
    labels = [
        f"{assay}-{condition}-{replicate}"
        for assay in ("RNA", "RIBO")
        for condition in ("A", "B")
        for replicate in (1, 2)
    ]
    feature_id = "NC_000913.3:200-499:+"
    (readcounts / "differential_expression_read_counts.csv").write_text(
        "Identifier," + ",".join(labels) + "\n" + feature_id + "," + ",".join(["20"] * 8) + "\n"
    )
    (readcounts / "independant_annotation.gff").write_text(
        "#hribo-gff-read-counts-v1\t...\n"
        "NC_000913.3\tHRIBO\tCDS\t200\t499\t.\t+\t.\tID=gene1;Name=gene1\t20\n"
    )
    for tool, header, values in (
        ("xtail", "log2FC_TE_final,pvalue_adjusted", "1.3,0.02"),
        ("riborex", "log2FC,pvalue_adjusted", "1.2,0.03"),
        (
            "deltate",
            "RNA_log2FC,RNA_pvalue_adjusted,RIBO_log2FC,RIBO_pvalue_adjusted,TE_log2FC,TE_pvalue_adjusted",
            "2,0.01,2.5,0.01,1.2,0.02",
        ),
    ):
        folder = workdir / tool
        folder.mkdir()
        (folder / f"{tool}_all.csv").write_text(
            f"gene_id,contrast,{header}\n{feature_id},{tool}_A-B,{values}\n"
        )

    result = subprocess.run(
        [
            *snakemake_command,
            "diffex_summary/condition_overview.html", "--cores", "1",
            "--allowed-rules", "conditionOverview",
            "--snakefile", str(REPO / "workflow" / "Snakefile"),
            "--directory", str(workdir), "--configfile", str(config_path),
        ],
        capture_output=True, text=True,
        env={**os.environ, "XDG_CACHE_HOME": str(workdir / ".cache")},
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (workdir / "diffex_summary" / "condition_matrix.tsv").is_file()
    assert (workdir / "diffex_summary" / "contrast_matrix.tsv").is_file()
    assert (workdir / "diffex_summary" / "browser_tracks.tsv").is_file()
    assert (workdir / "diffex_summary" / "browser" / "tracks_manifest.json").is_file()
