"""Runtime regression tests for staging FASTQs and compressed references."""

import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parent.parent
STAGER = REPO / "workflow" / "scripts" / "stage_input.py"
TRIMMING_RULES = (REPO / "workflow" / "rules" / "trimming.smk").read_text()
PREPROCESSING_RULES = (
    REPO / "workflow" / "rules" / "preprocessing.smk"
).read_text()


def run_stager(operation, source, destination, cwd=None):
    return subprocess.run(
        [sys.executable, str(STAGER), operation, str(source), str(destination)],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def test_relative_fastq_becomes_an_absolute_symlink(tmp_path):
    project = tmp_path / "project with spaces"
    source = project / "input reads" / "sample.fastq.gz"
    destination = project / "trimlink" / "sample.fastq.gz"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"fastq fixture")

    result = run_stager(
        "link",
        source.relative_to(project),
        destination.relative_to(project),
        cwd=project,
    )

    assert result.returncode == 0, result.stderr
    assert destination.is_symlink()
    assert Path(os.readlink(destination)) == source.absolute()
    assert destination.read_bytes() == source.read_bytes()


def test_absolute_fastq_is_not_prefixed_with_the_working_directory(tmp_path):
    project = tmp_path / "project"
    source = tmp_path / "external reads" / "sample.fastq.gz"
    destination = project / "trimlink" / "sample.fastq.gz"
    project.mkdir()
    source.parent.mkdir()
    source.write_bytes(b"fastq fixture")

    result = run_stager("link", source, destination.relative_to(project), cwd=project)

    assert result.returncode == 0, result.stderr
    assert Path(os.readlink(destination)) == source
    assert destination.read_bytes() == source.read_bytes()


def test_fastq_link_rerun_atomically_replaces_the_previous_target(tmp_path):
    first = tmp_path / "first.fastq.gz"
    second = tmp_path / "second.fastq.gz"
    destination = tmp_path / "trimlink" / "sample.fastq.gz"
    first.write_bytes(b"first")
    second.write_bytes(b"second")

    assert run_stager("link", first, destination).returncode == 0
    result = run_stager("link", second, destination)

    assert result.returncode == 0, result.stderr
    assert destination.is_symlink()
    assert Path(os.readlink(destination)) == second
    assert destination.read_bytes() == b"second"


def test_fastq_link_replaces_a_dangling_previous_target(tmp_path):
    source = tmp_path / "current.fastq.gz"
    destination = tmp_path / "trimlink" / "sample.fastq.gz"
    source.write_bytes(b"current")
    destination.parent.mkdir()
    destination.symlink_to(tmp_path / "deleted.fastq.gz")

    result = run_stager("link", source, destination)

    assert result.returncode == 0, result.stderr
    assert Path(os.readlink(destination)) == source
    assert destination.read_bytes() == b"current"


def test_fastq_link_rejects_a_missing_source_without_replacing_output(tmp_path):
    missing = tmp_path / "missing.fastq.gz"
    destination = tmp_path / "trimlink" / "sample.fastq.gz"
    destination.parent.mkdir()
    destination.write_bytes(b"keep me")

    result = run_stager("link", missing, destination)

    assert result.returncode != 0
    assert "does not exist or is not a file" in result.stderr
    assert destination.read_bytes() == b"keep me"


@pytest.mark.parametrize(
    "name, contents",
    [
        ("reference genome.fa.gz", b">chr1 description\nACGTACGT\n"),
        (
            "reference annotation.gff.gz",
            b"##gff-version 3\nchr1\ttest\tCDS\t1\t6\t.\t+\t0\tID=cds1\n",
        ),
    ],
)
def test_gzip_references_are_materialized_as_plain_text(tmp_path, name, contents):
    source = tmp_path / "compressed references" / name
    destination = tmp_path / "workflow outputs" / name.removesuffix(".gz")
    source.parent.mkdir()
    with gzip.open(source, "wb") as handle:
        handle.write(contents)

    result = run_stager("text", source, destination)

    assert result.returncode == 0, result.stderr
    assert destination.read_bytes() == contents
    assert destination.read_bytes()[:2] != b"\x1f\x8b"


def test_plain_reference_is_copied_unchanged(tmp_path):
    source = tmp_path / "reference.fa"
    destination = tmp_path / "outputs" / "genome.fa"
    contents = b">chr1\nACGT\n"
    source.write_bytes(contents)

    result = run_stager("text", source, destination)

    assert result.returncode == 0, result.stderr
    assert destination.read_bytes() == contents


def test_broken_gzip_does_not_replace_a_previous_reference(tmp_path):
    source = tmp_path / "broken.fa.gz"
    destination = tmp_path / "outputs" / "genome.fa"
    source.write_bytes(b"\x1f\x8bnot-a-complete-gzip-stream")
    destination.parent.mkdir()
    destination.write_bytes(b">previous\nACGT\n")

    result = run_stager("text", source, destination)

    assert result.returncode != 0
    assert destination.read_bytes() == b">previous\nACGT\n"
    assert list(destination.parent.glob(f".{destination.name}.*")) == []


def test_rules_delegate_to_the_quoted_stager():
    assert "os.getcwd() +" not in TRIMMING_RULES
    for placeholder in (
        "{input.stager:q}",
        "{input.fastq:q}",
        "{input.fastq1:q}",
        "{input.fastq2:q}",
    ):
        assert placeholder in TRIMMING_RULES

    assert "cp {input" not in PREPROCESSING_RULES
    assert PREPROCESSING_RULES.count("{input.stager:q}") == 2
    assert "{input.genome:q}" in PREPROCESSING_RULES
    assert "{input.annotation:q}" in PREPROCESSING_RULES
    assert "python3 {input.converter:q}" in PREPROCESSING_RULES
    assert "-o {output.annotation:q}" in PREPROCESSING_RULES
