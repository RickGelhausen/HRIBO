"""Executed workflow smoke tests using real rule-specific bioinformatics tools."""

from __future__ import annotations

import gzip
import os
import random
import subprocess
import sys
from pathlib import Path

import pysam
import yaml


REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow/Snakefile"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def overlapping_occurrences(sequence: str, query: str) -> int:
    return sum(
        sequence.startswith(query, start)
        for start in range(len(sequence) - len(query) + 1)
    )


def assert_uniquely_mappable(genome: str, query: str) -> None:
    assert (
        overlapping_occurrences(genome, query)
        + overlapping_occurrences(genome, reverse_complement(query))
        == 1
    )


def write_mapping_fixture(workdir: Path) -> tuple[Path, dict[str, tuple[int, bool]]]:
    inputs = workdir / "input data"
    inputs.mkdir(parents=True)

    rng = random.Random(20260902)
    sequence = "".join(rng.choice("ACGT") for _ in range(6000))
    planted = {
        "retained-forward": (sequence[1000:1030], 1000, False),
        "retained-reverse": (reverse_complement(sequence[2000:2030]), 2000, True),
        "removed-rrna": (sequence[300:330], 300, False),
        "removed-trna": (reverse_complement(sequence[5000:5030]), 5000, True),
    }
    for query, _, _ in planted.values():
        assert_uniquely_mappable(sequence, query)

    genome_contents = f">fixture_contig synthetic mapping fixture\n{sequence}\n"
    genome = inputs / "tiny genome.fa.gz"
    with gzip.open(genome, "wt") as handle:
        handle.write(genome_contents)

    annotation_contents = (
        "##gff-version 3\n"
        "fixture_contig\tfixture\trRNA\t201\t500\t.\t+\t.\tID=rrna1\n"
        "fixture_contig\tfixture\tCDS\t801\t2498\t.\t+\t0\t"
        "ID=cds1;locus_tag=fixture_cds\n"
        "fixture_contig\tfixture\ttRNA\t4901\t5200\t.\t-\t.\tID=trna1\n"
    )
    annotation = inputs / "tiny annotation.gff.gz"
    with gzip.open(annotation, "wt") as handle:
        handle.write(annotation_contents)

    fastq = inputs / "reads.fastq.gz"
    with gzip.open(fastq, "wt") as handle:
        for name, (query, _, _) in planted.items():
            handle.write(f"@{name}\n{query}\n+\n{'I' * len(query)}\n")

    samples = inputs / "samples.tsv"
    samples.write_text(
        "method\tcondition\treplicate\tfastqFile\tfastqFile2\n"
        "RIBO\tSmoke\t1\tinput data/reads.fastq.gz\t\n"
    )

    config = yaml.safe_load((REPO / "config/config.yaml").read_text())
    config["biologySettings"].update(
        {
            "genome": str(genome),
            "annotation": str(annotation),
            "samples": str(samples),
        }
    )
    config["predictionSettings"]["deepribo"] = "off"
    config["workflowSettings"]["stages"] = ["mapping"]
    config_path = inputs / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    truth = {
        name: (position, reverse)
        for name, (_, position, reverse) in planted.items()
    }
    return config_path, truth


def test_mapping_stage_executes_with_real_tools_and_filters_structural_rna(
    snakemake_command,
    tmp_path,
):
    workdir = tmp_path / "mapping run with spaces"
    workdir.mkdir()
    config_path, truth = write_mapping_fixture(workdir)

    configured_prefix = os.environ.get("HRIBO_TEST_CONDA_PREFIX")
    conda_prefix = (
        Path(configured_prefix)
        if configured_prefix
        else REPO / ".snakemake/conda"
    )
    conda_prefix = conda_prefix.resolve()

    environment = {
        **os.environ,
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
            "--show-failed-logs",
            "--rerun-incomplete",
            "--notemp",
            "--software-deployment-method",
            "conda",
            "--conda-prefix",
            str(conda_prefix),
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
        timeout=1200,
    )
    rendered = result.stdout + result.stderr
    assert result.returncode == 0, rendered

    staged_genome = workdir / "genomes/genome.fa"
    staged_annotation = workdir / "annotation/annotation.gff"
    assert staged_genome.read_bytes()[:2] != b"\x1f\x8b"
    assert staged_annotation.read_bytes()[:2] != b"\x1f\x8b"
    assert "fixture_contig" in staged_genome.read_text()
    assert "ID=rrna1" in staged_annotation.read_text()
    assert (workdir / "annotation/rrna.bed").read_text().splitlines() == [
        "fixture_contig\t200\t500\t.\t.\t+",
        "fixture_contig\t4900\t5200\t.\t.\t-",
    ]

    prefilter_bam = workdir / "rRNAbam/RIBO-Smoke-1.bam"
    with pysam.AlignmentFile(prefilter_bam, "rb") as alignment:
        prefilter_records = list(alignment.fetch(until_eof=True))
    assert len(prefilter_records) == 4
    assert {record.query_name for record in prefilter_records} == set(truth)
    for record in prefilter_records:
        expected_position, expected_reverse = truth[record.query_name]
        assert record.reference_start == expected_position
        assert record.is_reverse is expected_reverse
        assert record.get_tag("NH") == 1

    bam = workdir / "maplink/RIBO-Smoke-1.bam"
    index = Path(f"{bam}.bai")
    assert bam.is_symlink()
    assert not os.path.isabs(os.readlink(bam))
    assert index.is_file()

    with pysam.AlignmentFile(bam, "rb") as alignment:
        records = list(alignment.fetch("fixture_contig"))
    assert len(records) == 2
    assert {record.query_name for record in records} == {
        "retained-forward",
        "retained-reverse",
    }
    for record in records:
        expected_position, expected_reverse = truth[record.query_name]
        assert record.reference_start == expected_position
        assert record.is_reverse is expected_reverse
        assert record.get_tag("NH") == 1

    moved = tmp_path / "moved mapping result"
    workdir.rename(moved)
    moved_bam = moved / "maplink/RIBO-Smoke-1.bam"
    with pysam.AlignmentFile(moved_bam, "rb") as alignment:
        moved_records = list(alignment.fetch("fixture_contig"))
        assert len(moved_records) == 2
        assert {record.query_name for record in moved_records} == {
            "retained-forward",
            "retained-reverse",
        }
