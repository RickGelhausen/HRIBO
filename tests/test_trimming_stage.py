"""Executed and dry-run contracts for the standalone trimming stage."""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
from pathlib import Path

import yaml


REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"


def _write_fastq(path: Path, name: str, sequence: str) -> str:
    contents = f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"
    with gzip.open(path, "wt") as handle:
        handle.write(contents)
    return contents


def _trimming_project(tmp_path: Path) -> tuple[Path, Path, dict[str, str]]:
    workdir = tmp_path / "trimming workflow with spaces"
    inputs = workdir / "input reads"
    inputs.mkdir(parents=True)

    single = inputs / "single reads.fastq.gz"
    paired_forward = inputs / "paired forward.fastq.gz"
    paired_reverse = inputs / "paired reverse.fastq.gz"
    expected = {
        "single": _write_fastq(single, "single", "ACGTACGTACGT"),
        "paired_forward": _write_fastq(
            paired_forward, "paired-forward", "AAAACCCCGGGG"
        ),
        "paired_reverse": _write_fastq(
            paired_reverse, "paired-reverse", "TTTTGGGGCCCC"
        ),
    }

    samples = inputs / "samples.tsv"
    samples.write_text(
        "method\tcondition\treplicate\tfastqFile\tfastqFile2\n"
        f"RNA\tSingle\t1\t{single}\t\n"
        f"RNA\tPaired\t1\t{paired_forward}\t{paired_reverse}\n"
    )

    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["biologySettings"].update(
        {
            "genome": str(inputs / "intentionally absent genome.fa"),
            "annotation": str(inputs / "intentionally absent annotation.gff"),
            "samples": str(samples),
        }
    )
    config["predictionSettings"]["deepribo"] = "off"
    config["workflowSettings"]["stages"] = ["trimming"]
    config_path = inputs / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return workdir, config_path, expected


def _environment(workdir: Path, tools: Path | None = None) -> dict[str, str]:
    path_entries = [str(Path(sys.executable).parent)]
    if tools is not None:
        path_entries.insert(0, str(tools))
    path_entries.append(os.environ.get("PATH", ""))
    return {
        **os.environ,
        "PATH": os.pathsep.join(path_entries),
        "XDG_CACHE_HOME": str(workdir / ".cache"),
    }


def _write_executable(path: Path, source: str) -> None:
    path.write_text(source)
    path.chmod(0o755)


def _fake_trimming_tools(workdir: Path) -> Path:
    tools = workdir / "fake tools"
    tools.mkdir()

    _write_executable(
        tools / "cutadapt",
        """#!/usr/bin/env python3
import gzip
import sys
from pathlib import Path

arguments = sys.argv[1:]
outputs = [arguments[arguments.index("-o") + 1]]
if "-p" in arguments:
    outputs.append(arguments[arguments.index("-p") + 1])
sources = arguments[-len(outputs):]
for source, output in zip(sources, outputs):
    with gzip.open(source, "rt") as handle:
        contents = handle.read()
    Path(output).write_text(contents)
""",
    )
    _write_executable(
        tools / "pear",
        """#!/usr/bin/env python3
import sys
from pathlib import Path

arguments = sys.argv[1:]
forward = Path(arguments[arguments.index("-f") + 1]).read_text()
reverse = Path(arguments[arguments.index("-r") + 1]).read_text()
prefix = arguments[arguments.index("-o") + 1]
Path(prefix + ".assembled.fastq").write_text(forward + reverse)
""",
    )
    _write_executable(
        tools / "fastqc",
        """#!/usr/bin/env python3
import sys
from pathlib import Path

arguments = sys.argv[1:]
outdir = Path(arguments[arguments.index("-o") + 1])
name = Path(arguments[-1]).name
for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
    if name.endswith(suffix):
        name = name[:-len(suffix)]
        break
(outdir / f"{name}_fastqc.html").write_text("fake FastQC HTML\\n")
(outdir / f"{name}_fastqc.zip").write_text("fake FastQC archive\\n")
""",
    )
    return tools


def test_trimming_stage_dag_publishes_single_and_paired_processed_reads(
    snakemake_command,
    tmp_path,
):
    workdir, config_path, _ = _trimming_project(tmp_path)
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
            str(workdir),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env=_environment(workdir),
    )
    rendered = result.stdout + result.stderr

    assert result.returncode == 0, rendered
    for rule in (
        "trim_single",
        "trim_paired",
        "merge_fastq",
        "fastqcraw_single",
        "fastqcraw_paired",
        "fastqctrimmed_single",
        "fastqctrimmed_paired",
    ):
        assert f"rule {rule}" in rendered
    assert "trimmed/RNA-Single-1.fastq" in rendered
    assert "trimmed/RNA-Paired-1.fastq" in rendered
    assert "rule retrieveGenome" not in rendered
    assert "rule checkAnnotation" not in rendered
    assert "rule map:" not in rendered


def test_trimming_stage_execution_retains_mapping_ready_reads_for_both_layouts(
    snakemake_command,
    tmp_path,
):
    workdir, config_path, expected = _trimming_project(tmp_path)
    tools = _fake_trimming_tools(workdir)
    result = subprocess.run(
        [
            *snakemake_command,
            "all",
            "--cores",
            "1",
            "--printshellcmds",
            "--show-failed-logs",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env=_environment(workdir, tools),
        timeout=300,
    )
    rendered = result.stdout + result.stderr

    assert result.returncode == 0, rendered
    single = workdir / "trimmed" / "RNA-Single-1.fastq"
    paired = workdir / "trimmed" / "RNA-Paired-1.fastq"
    assert single.read_text() == expected["single"]
    assert paired.read_text() == (
        expected["paired_forward"] + expected["paired_reverse"]
    )

    expected_reports = {
        "qc/1raw/RNA-Single-1-raw_fastqc.html",
        "qc/1raw/RNA-Paired-1-raw-q_fastqc.html",
        "qc/1raw/RNA-Paired-1-raw-p_fastqc.html",
        "qc/2trimmed/RNA-Single-1-trimmed_fastqc.html",
        "qc/2trimmed/RNA-Paired-1-trimmed_q_fastqc.html",
        "qc/2trimmed/RNA-Paired-1-trimmed_p_fastqc.html",
    }
    assert all((workdir / report).is_file() for report in expected_reports)
    assert not (workdir / "genomes").exists()
    assert not (workdir / "annotation").exists()
