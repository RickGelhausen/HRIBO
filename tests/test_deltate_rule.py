"""Regression tests for deltaTE's Snakemake and runtime failure contract."""

import os
import stat
import subprocess
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parent.parent
RULE = (REPO / "workflow" / "rules" / "diffex_deltate.smk").read_text()
RUNNER = REPO / "workflow" / "scripts" / "run_deltate.sh"


@pytest.fixture
def deltate_run(tmp_path):
    """Create input/output paths plus a controllable fake ``DTEG.R``."""
    bin_dir = tmp_path / "fake bin"
    bin_dir.mkdir()
    executable = bin_dir / "DTEG.R"
    executable.write_text(
        """#!/usr/bin/env bash
set -euo pipefail
mode=${FAKE_DELTATE_MODE:-success}
result_dir=${5%/}
mkdir -p "$result_dir/fold_changes"

if [[ "$mode" == "fail" ]]; then
    echo "simulated DTEG failure" >&2
    exit 23
fi

printf 'baseMean\\tlog2FoldChange\\tlfcSE\\tpvalue\\tpadj\\n' > "$result_dir/fold_changes/deltaRibo.txt"
printf 'baseMean\\tlog2FoldChange\\tlfcSE\\tpvalue\\tpadj\\n' > "$result_dir/fold_changes/deltaRNA.txt"
if [[ "$mode" != "missing" ]]; then
    printf 'baseMean\\tlog2FoldChange\\tlfcSE\\tstat\\tpvalue\\tpadj\\n' > "$result_dir/fold_changes/deltaTE.txt"
fi
if [[ "$mode" != "header_only" ]]; then
    printf 'gene_1\\t10\\t1\\t0.2\\t0.01\\t0.02\\n' >> "$result_dir/fold_changes/deltaRibo.txt"
    printf 'gene_1\\t10\\t0.5\\t0.2\\t0.01\\t0.02\\n' >> "$result_dir/fold_changes/deltaRNA.txt"
    if [[ "$mode" != "missing" ]]; then
        printf 'gene_1\\t10\\t0.5\\t0.2\\t2.5\\t0.01\\t0.02\\n' >> "$result_dir/fold_changes/deltaTE.txt"
    fi
fi
if [[ "$mode" == "bad_pdf" ]]; then
    printf 'not a PDF\\n' > "$result_dir/Result_figures.pdf"
else
    printf '%%PDF-1.4\\n' > "$result_dir/Result_figures.pdf"
fi
"""
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)

    inputs = []
    for name in ("ribo counts.txt", "rna counts.txt", "samples info.txt"):
        path = tmp_path / name
        path.write_text("fixture\n")
        inputs.append(path)

    result_dir = tmp_path / "result directory"
    outputs = [
        result_dir / "fold_changes" / "deltaRibo.txt",
        result_dir / "fold_changes" / "deltaRNA.txt",
        result_dir / "fold_changes" / "deltaTE.txt",
    ]
    figure_source = result_dir / "Result_figures.pdf"
    figure_output = tmp_path / "published figures.pdf"
    command = [
        "bash",
        str(RUNNER),
        *(str(path) for path in inputs),
        str(result_dir),
        *(str(path) for path in outputs),
        str(figure_source),
        str(figure_output),
        str(executable),
    ]
    environment = os.environ.copy()

    return command, environment, outputs, figure_source, figure_output


def execute(fixture, mode):
    command, environment, *_ = fixture
    environment = {**environment, "FAKE_DELTATE_MODE": mode}
    return subprocess.run(command, capture_output=True, text=True, env=environment)


def test_deltate_rule_delegates_to_the_checked_runner():
    assert 'runner=str(SCRIPTS / "run_deltate.sh")' in RULE
    assert 'patcher=str(SCRIPTS / "patch_deltate.py")' in RULE
    assert "rule prepareDeltaTEScript:" in RULE
    assert "engine=rules.prepareDeltaTEScript.output.script" in RULE
    assert "{input.engine:q}" in RULE
    assert "bash {input.runner:q}" in RULE
    assert "|| true" not in RULE
    assert "touch {output" not in RULE
    for output in ("deltaRibo.txt", "deltaRNA.txt", "deltaTE.txt", "_figures.pdf"):
        declaration = next(line for line in RULE.splitlines() if output in line)
        assert "ensure(" in declaration
        assert "non_empty=True" in declaration


def test_deltate_propagates_the_tool_exit_and_removes_stale_outputs(deltate_run):
    _, _, outputs, figure_source, figure_output = deltate_run
    for path in (*outputs, figure_source, figure_output):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("stale\n")

    result = execute(deltate_run, "fail")

    assert result.returncode == 23
    assert "simulated DTEG failure" in result.stderr
    assert not any(path.exists() for path in (*outputs, figure_source, figure_output))


def test_deltate_rejects_a_success_exit_with_missing_outputs(deltate_run):
    result = execute(deltate_run, "missing")

    assert result.returncode != 0
    assert "completed without required non-empty output" in result.stderr


def test_deltate_validates_tables_and_pdf_before_publication(deltate_run):
    header_only = execute(deltate_run, "header_only")
    assert header_only.returncode != 0
    assert "has no valid data rows" in header_only.stderr

    bad_pdf = execute(deltate_run, "bad_pdf")
    assert bad_pdf.returncode != 0
    assert "figure is not a PDF" in bad_pdf.stderr

    result = execute(deltate_run, "success")
    _, _, outputs, figure_source, figure_output = deltate_run

    assert result.returncode == 0, result.stderr
    assert all(path.stat().st_size > 0 for path in outputs)
    assert figure_output.read_bytes() == figure_source.read_bytes()
