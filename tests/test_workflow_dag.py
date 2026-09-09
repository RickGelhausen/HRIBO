"""Fixture-backed DAG construction checks for the supported workflow breadth."""

from __future__ import annotations

import gzip
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml


REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow/Snakefile"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def dry_run(
    *,
    stages,
    snakemake_command,
    genome_file,
    annotation_file,
    samples,
    tmp_path,
    deepribo=None,
):
    sample_path = tmp_path / "samples.tsv"
    samples.fillna("").to_csv(sample_path, sep="\t", index=False)

    workflow_config = yaml.safe_load((REPO / "config/config.yaml").read_text())
    workflow_config["biologySettings"].update(
        {
            "genome": str(genome_file),
            "annotation": str(annotation_file),
            "samples": str(sample_path),
        }
    )
    workflow_config["workflowSettings"]["stages"] = stages
    if deepribo is not None:
        workflow_config["predictionSettings"]["deepribo"] = deepribo

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
            str(SNAKEFILE),
            "--directory",
            str(tmp_path),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / ".cache")},
    )
    rendered = result.stdout + result.stderr
    assert result.returncode == 0, rendered
    return rendered


def test_mapping_stage_constructs_complete_fixture_dag(
    snakemake_command,
    genome_file,
    annotation_file,
    samples,
    tmp_path,
):
    rendered = dry_run(
        stages=["mapping"],
        snakemake_command=snakemake_command,
        genome_file=genome_file,
        annotation_file=annotation_file,
        samples=samples,
        tmp_path=tmp_path,
    )

    assert "rule genomeSegemehlIndex" in rendered
    assert "rule map" in rendered
    assert "rule maplink" in rendered
    assert "maplink/RIBO-A-1.bam" in rendered
    assert "maplink/RNA-B-2.bam.bai" in rendered
    assert "rule createOverviewTable" not in rendered
    assert "rule updatedAnnotation" not in rendered
    assert "tracks/updated_annotation.gff" not in rendered


def test_full_preset_constructs_every_analysis_branch(
    snakemake_command,
    genome_file,
    annotation_file,
    samples,
    tmp_path,
):
    alphanumeric_samples = samples.copy()
    alphanumeric_samples["condition"] = alphanumeric_samples["condition"].replace(
        {"A": "A1"}
    )
    rendered = dry_run(
        stages="full",
        snakemake_command=snakemake_command,
        genome_file=genome_file,
        annotation_file=annotation_file,
        samples=alphanumeric_samples,
        tmp_path=tmp_path,
    )

    for rule in (
        "multiqc",
        "metageneProfiling",
        "tisAdvisor",
        "predictDeepRibo",
        "prepareReparationAnnotation",
        "reparation",
        "prepareDeltaTEScript",
        "xtail",
        "riborex",
        "deltate",
        "createOverviewTable",
        "updatedAnnotation",
    ):
        assert f"rule {rule}" in rendered
    assert "logs/A1-1_reparation.log" in rendered
    assert "reparation/annotation.gtf" in rendered
    assert "deltate/DTEG.R" in rendered
    assert "tracks/updated_annotation.gff" in rendered
    assert "build_updated_annotation.py" in rendered
    assert "tracks/deepribo_merged_plus.gff" in rendered
    for overview_output in (
        "auxiliary/overview.xlsx",
        "auxiliary/overview.tsv",
        "auxiliary/overview.gff",
        "auxiliary/overview_misc.gff",
    ):
        assert overview_output in rendered

    # Container jobs run with the selected project directory as their home.
    # Bundled helpers therefore have to resolve through Snakemake's source cache,
    # which it bind-mounts independently of where the HRIBO checkout lives.
    for helper in (
        "patch_deltate.py",
        "run_deltate.sh",
        "run_reparation.py",
        "patch_deepribo_scurve.py",
        "deepribo_data_parser.py",
        "run_parameter_estimation.py",
        "parameter_estimation.R",
    ):
        assert f"{helper} (cached)" in rendered
        assert re.search(rf"source-cache/\S*{re.escape(helper)}", rendered)


@pytest.mark.parametrize(
    ("deepribo", "uses_deepribo"),
    [("on", True), ("off", False)],
    ids=["deepribo-on", "deepribo-off"],
)
def test_predictions_stage_builds_the_updated_annotation_in_both_modes(
    deepribo,
    uses_deepribo,
    snakemake_command,
    genome_file,
    annotation_file,
    samples,
    tmp_path,
):
    rendered = dry_run(
        stages=["predictions"],
        snakemake_command=snakemake_command,
        genome_file=genome_file,
        annotation_file=annotation_file,
        samples=samples,
        tmp_path=tmp_path,
        deepribo=deepribo,
    )

    assert "rule updatedAnnotation" in rendered
    assert "tracks/updated_annotation.gff" in rendered
    assert "tracks/reparation_annotated.gff" in rendered
    assert ("tracks/deepribo_merged_plus.gff" in rendered) is uses_deepribo
    assert ("rule predictDeepRibo" in rendered) is uses_deepribo


def test_obsolete_updated_annotation_chain_is_absent():
    production_files = [SNAKEFILE, *(REPO / "workflow" / "rules").glob("*.smk")]
    production_text = "\n".join(path.read_text() for path in production_files)

    for obsolete in (
        "tracks/totalAnnotation.gff",
        "uniteAnnotation",
        "newAnnotationDeepRibo",
        "newAnnotationReparationOnly",
        "annotation_unite.py",
    ):
        assert obsolete not in production_text
    assert not (REPO / "workflow" / "rules" / "conditionals.smk").exists()
    assert not (REPO / "workflow" / "scripts" / "annotation_unite.py").exists()


def test_genome_track_stage_executes_real_workflow_with_compressed_spaced_input(
    snakemake_command,
    tmp_path,
):
    workdir = tmp_path / "executed workflow with spaces"
    inputs = workdir / "input data"
    inputs.mkdir(parents=True)

    # Each requested motif occurs once on each strand. Ns prevent accidental
    # cross-boundary motifs and remain valid IUPAC sequence characters.
    sequence = (
        "ATGNNNGTGNNNTTGNNNTAGNNNTGANNNTAANNNAAGGNNN"
        "CATNNNCACNNNCAANNNCTANNNTCANNNTTANNNCCTT"
    )
    genome = inputs / "tiny genome.fa.gz"
    with gzip.open(genome, "wt") as handle:
        handle.write(f">tiny_contig integration fixture\n{sequence}\n")

    samples = inputs / "samples.tsv"
    samples.write_text(
        "method\tcondition\treplicate\tfastqFile\tfastqFile2\n"
        "RNA\tA1\t1\tnot-needed.fastq.gz\t\n"
    )

    workflow_config = yaml.safe_load((REPO / "config/config.yaml").read_text())
    workflow_config["biologySettings"].update(
        {
            "genome": str(genome),
            "annotation": str(inputs / "not-needed.gff"),
            "samples": str(samples),
            "alternativeStartCodons": ["GTG", "TTG"],
        }
    )
    workflow_config["predictionSettings"]["deepribo"] = "off"
    workflow_config["workflowSettings"]["stages"] = ["genome_tracks"]
    config_path = inputs / "config.yaml"
    config_path.write_text(yaml.safe_dump(workflow_config, sort_keys=False))

    environment = {
        **os.environ,
        # Rule commands use python3. Keep them on the same complete Python
        # environment that is running pytest, even when Snakemake itself is
        # reached through `conda run` in a developer checkout.
        "PATH": os.pathsep.join(
            [str(Path(sys.executable).parent), os.environ.get("PATH", "")]
        ),
        "XDG_CACHE_HOME": str(workdir / ".cache"),
    }
    result = subprocess.run(
        [
            *snakemake_command,
            "all",
            "--cores",
            "1",
            "--printshellcmds",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env=environment,
        timeout=300,
    )
    rendered = result.stdout + result.stderr
    assert result.returncode == 0, rendered

    assert (workdir / "genomes/genome.fa").read_text() == (
        f">tiny_contig integration fixture\n{sequence}\n"
    )
    expected_names = {
        "potentialStartCodons.gff": {"ATG"},
        "potentialAlternativeStartCodons.gff": {"GTG", "TTG"},
        "potentialStopCodons.gff": {"TAG", "TGA", "TAA"},
        "potentialRibosomeBindingSite.gff": {"AAGG"},
    }
    for filename, motifs in expected_names.items():
        lines = [
            line.split("\t")
            for line in (workdir / "tracks" / filename).read_text().splitlines()
            if line and not line.startswith("#")
        ]
        assert len(lines) == len(motifs) * 2
        assert {line[0] for line in lines} == {"tiny_contig"}
        assert {line[6] for line in lines} == {"+", "-"}
        names = {
            field.split("=", 1)[1]
            for line in lines
            for field in line[8].rstrip(";").split(";")
            if field.startswith("Name=")
        }
        assert names == motifs
        for line in lines:
            name = next(
                field.split("=", 1)[1]
                for field in line[8].rstrip(";").split(";")
                if field.startswith("Name=")
            )
            genomic_segment = sequence[int(line[3]) - 1 : int(line[4])]
            observed_motif = (
                genomic_segment
                if line[6] == "+"
                else reverse_complement(genomic_segment)
            )
            assert observed_motif == name
