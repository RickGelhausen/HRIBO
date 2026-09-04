"""Regression coverage for DeepRibo S-curve estimation boundaries."""

import base64
import hashlib
import os
import signal
import shutil
import stat
import subprocess
import sys
import time
import zlib
from pathlib import Path

import pytest

from patch_deepribo_scurve import (
    EXPECTED_SHA256,
    PatchError,
    materialize_patched_script,
)


REPO = Path(__file__).resolve().parent.parent
WRAPPER = REPO / "workflow" / "scripts" / "parameter_estimation.R"
LAUNCHER = REPO / "workflow" / "scripts" / "run_parameter_estimation.py"
RULE = (REPO / "workflow" / "rules" / "deepribo.smk").read_text()

SOURCE_FRAGMENT = """#!/usr/bin/env Rscript
get_cutoff_values <- function(path,dest){
  while (MINCOV>0.60){
    MINCOV <- round(predict(bent_curve, MINRPKM), digits = 6)
  }
  plot(df$rpk_elo, df$coverage_elo, xlab="RPKM",ylab="Coverage")
}
"""

VALID_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk"
    "+A8AAQUBAScY42YAAAAASUVORK5CYII="
)


def add_png_text_chunk(png, text):
    """Return a distinct valid PNG without needing an image dependency."""
    chunk_type = b"tEXt"
    payload = text.encode()
    chunk = (
        len(payload).to_bytes(4, "big")
        + chunk_type
        + payload
        + (zlib.crc32(chunk_type + payload) & 0xFFFFFFFF).to_bytes(4, "big")
    )
    return png[:-12] + chunk + png[-12:]


CURRENT_PNG = add_png_text_chunk(VALID_PNG, "publication=current")


def write_curve(path, rows=5):
    values = ["label,rpk_elo,coverage_elo"]
    for index in range(rows):
        values.append(f"1,{index + 1},{0.1 + index * 0.1}")
    path.write_text("\n".join(values) + "\n")


def wrapper_command(data, output, plot, engine, receipt=None):
    command = [
        sys.executable,
        str(LAUNCHER),
        str(WRAPPER),
        "--file",
        str(data),
        "--out",
        str(output),
        "--plot",
        str(plot),
    ]
    if receipt is not None:
        command.extend(["--receipt", str(receipt)])
    command.extend(["--engine", str(engine)])
    return command


def run_wrapper(
    data, output, plot, engine, environment=None, cwd=None, receipt=None
):
    return subprocess.run(
        wrapper_command(data, output, plot, engine, receipt=receipt),
        capture_output=True,
        text=True,
        env=environment,
        cwd=cwd,
    )


def transaction_path(output, state):
    return output.parent / f".{output.name}.pair-transaction.{state}"


def transaction_backup(target, index):
    return target.parent / f".{target.name}.pair-transaction-old-{index}"


def initialize_transaction(path, targets, existed=(True, True)):
    path.mkdir()
    for index, was_present in enumerate(existed, 1):
        state = "present" if was_present else "absent"
        (path / f"{state}-{index}").touch()
    for index, target in enumerate(targets, 1):
        (path / f"target-{index}").write_text(str(target.resolve()) + "\n")
    (path / "ready").touch()


def write_failing_engine(path):
    path.write_text(
        "get_cutoff_values <- function(path, dest) "
        "stop('failure after startup reconciliation')\n"
    )


def production_parameter_rule_run(tmp_path, snakemake_command):
    """Stage the included production rule around a controllable R engine."""

    workflow_root = tmp_path / "production parameter workflow"
    parsed = workflow_root / "deepribo/parsed/A-1"
    parsed.mkdir(parents=True)
    write_curve(parsed / "data_list.csv")
    engine = workflow_root / "deepribo/s_curve_cutoff_estimation.R"
    engine.write_text(
        r'''base_file_rename <- base::file.rename
assign("file.rename", function(from, to) {
  renamed <- base_file_rename(from, to)
  notification <- Sys.getenv("HRIBO_TEST_PARAMETER_INSTALL_NOTIFY")
  if (isTRUE(renamed) && nzchar(notification) && basename(to) == ".complete") {
    writeLines(as.character(Sys.getpid()), notification)
    while (TRUE) Sys.sleep(1)
  }
  renamed
}, envir=.GlobalEnv)

get_cutoff_values <- function(path, dest) {
  mode <- Sys.getenv("FAKE_PARAMETER_MODE", "first")
  if (mode == "fail") stop("simulated parameter-estimation failure")
  png(paste0(dest, ".png"), width=600, height=600)
  if (mode == "first") plot(1, 1) else plot(2, 1)
  invisible(dev.off())
  if (mode == "first") {
    list(min_RPKM=0.125, min_coverage=0.4)
  } else {
    list(min_RPKM=0.25, min_coverage=0.5)
  }
}
'''
    )
    for path, content in (
        (workflow_root / "genomes/genome.fa", ">chr1\nA\n"),
        (workflow_root / "annotation/annotation.gff", "##gff-version 3\n"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)

    snakefile = workflow_root / "ProductionParameterSnakefile"
    snakefile.write_text(
        f'''from pathlib import Path
import pandas as pd

SCRIPTS = Path({str(REPO / "workflow/scripts")!r})
conditions = ["A"]
samples = pd.DataFrame(
    [{{"method": "RIBO", "condition": "A", "replicate": "1"}}]
)
config["predictionSettings"] = {{"deepriboASiteOffset": 12}}


rule retrieveGenome:
    output:
        "genomes/genome.fa"


rule checkAnnotation:
    output:
        "annotation/annotation.gff"


include: {str(REPO / "workflow/rules/deepribo.smk")!r}
'''
    )
    receipt = workflow_root / "deepribo/cutoffs/A-1/.complete"
    invocation = [
        *snakemake_command,
        receipt.relative_to(workflow_root).as_posix(),
        "--snakefile",
        str(snakefile),
        "--cores",
        "1",
        "--allowed-rules",
        "parameterEstimation",
        "--rerun-incomplete",
        "--printshellcmds",
        "--show-failed-logs",
    ]
    environment = {
        **os.environ,
        "XDG_CACHE_HOME": str(workflow_root / ".cache"),
    }
    return invocation, environment, workflow_root, receipt, snakefile


def test_scurve_patcher_is_checksum_guarded_and_atomic(tmp_path):
    source = tmp_path / "installed estimator.R"
    output = tmp_path / "generated" / "estimator.R"
    source.write_text(SOURCE_FRAGMENT)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()

    materialize_patched_script(source, output, expected_sha256=digest)

    patched = output.read_text()
    assert "if (fit_idx < 5L)" in patched
    assert "non-finite bend" in patched
    assert "Log mean A-site occupancy per nucleotide" in patched
    assert "xlab=\"RPKM\"" not in patched
    assert output.stat().st_mode & stat.S_IXUSR

    output.write_text("previous valid estimator\n")
    with pytest.raises(PatchError, match="unexpected DeepRibo estimator"):
        materialize_patched_script(source, output, expected_sha256="0" * 64)
    assert output.read_text() == "previous valid estimator\n"

    drifted = SOURCE_FRAGMENT.replace("  while (MINCOV>0.60){", "  while (TRUE){")
    source.write_text(drifted)
    drifted_digest = hashlib.sha256(source.read_bytes()).hexdigest()
    with pytest.raises(PatchError, match="expected 1 occurrences"):
        materialize_patched_script(
            source, output, expected_sha256=drifted_digest
        )
    assert output.read_text() == "previous valid estimator\n"
    assert EXPECTED_SHA256 in (
        REPO / "workflow" / "scripts" / "patch_deepribo_scurve.py"
    ).read_text()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_validates_and_atomically_publishes_pair(tmp_path):
    directory = tmp_path / "paths with spaces"
    directory.mkdir()
    data = directory / "data list.csv"
    output = directory / "parameters.txt"
    plot = directory / "s curve.png"
    engine = directory / "fake estimator.R"
    write_curve(data)
    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1:5, 1:5)
  invisible(dev.off())
  list(min_RPKM=0.125, min_coverage=0.4)
}
"""
    )

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode == 0, result.stderr
    assert output.read_text() == "0.125,0.40000000000000002\n"
    assert plot.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert list(directory.glob(".parameters.txt.*")) == []
    assert list(directory.glob(".s curve.png.*")) == []


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_publishes_relative_output_paths(tmp_path):
    workdir = tmp_path / "work"
    output_directory = workdir / "deepribo" / "sample"
    output_directory.mkdir(parents=True)
    data = workdir / "data.csv"
    engine = workdir / "estimator.R"
    write_curve(data)
    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.125, min_coverage=0.4)
}
"""
    )
    output = Path("deepribo/sample/parameters.txt")
    plot = Path("deepribo/sample/s_curve.png")

    result = run_wrapper(
        Path("data.csv"), output, plot, Path("estimator.R"), cwd=workdir
    )

    assert result.returncode == 0, result.stderr
    assert (workdir / output).read_text() == "0.125,0.40000000000000002\n"
    assert (workdir / plot).read_bytes().startswith(b"\x89PNG")
    assert not transaction_path(workdir / output, "staging").exists()
    assert not transaction_path(workdir / output, "committed").exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_launcher_ignores_a_hostile_user_r_profile(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "estimator.R"
    profile = tmp_path / "hostile user profile.R"
    profile_marker = tmp_path / "profile-loaded"
    write_curve(data)
    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.125, min_coverage=0.4)
}
"""
    )
    profile.write_text(
        f'writeLines("loaded", {str(profile_marker)!r})\n'
        'stop("hostile R profile must not execute")\n'
    )

    result = run_wrapper(
        data,
        output,
        plot,
        engine,
        {**os.environ, "R_PROFILE_USER": str(profile)},
    )

    assert result.returncode == 0, result.stderr
    assert not profile_marker.exists()
    assert output.read_text() == "0.125,0.40000000000000002\n"
    assert plot.read_bytes().startswith(b"\x89PNG")


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_launcher_removes_inherited_r_package_overrides(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "estimator.R"
    hostile_library = tmp_path / "host supplied R library"
    hostile_library.mkdir()
    write_curve(data)
    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  hostile <- normalizePath(Sys.getenv("HRIBO_HOSTILE_R_LIBRARY"), mustWork=TRUE)
  if (hostile %in% .libPaths()) stop("host R library leaked into runtime")
  if (nzchar(Sys.getenv("R_ARCH"))) stop("host R architecture leaked into runtime")
  expected <- c("package:stats", "package:utils", "package:grDevices")
  if (!all(expected %in% search())) stop("default R packages were overridden")
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.125, min_coverage=0.4)
}
"""
    )
    environment = {
        **os.environ,
        "HRIBO_HOSTILE_R_LIBRARY": str(hostile_library),
        "R_ARCH": "/definitely-not-the-container-r-architecture",
        "R_DEFAULT_PACKAGES": "NULL",
        "R_HOME": "/definitely-not-the-container-r-home",
        "R_LIBS": str(hostile_library),
        "R_LIBS_SITE": str(hostile_library),
        "R_LIBS_USER": str(hostile_library),
    }

    result = run_wrapper(data, output, plot, engine, environment)

    assert result.returncode == 0, result.stderr
    assert "ignoring environment value of R_HOME" not in result.stderr
    assert output.read_text() == "0.125,0.40000000000000002\n"
    assert plot.read_bytes().startswith(b"\x89PNG")


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_rejects_sparse_input_before_loading_engine(tmp_path):
    data = tmp_path / "sparse.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "must not load.R"
    marker = tmp_path / "engine-loaded"
    write_curve(data, rows=4)
    engine.write_text(
        f'writeLines("loaded", {str(marker)!r})\n'
        "get_cutoff_values <- function(path, dest) stop('should not run')\n"
    )

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "requires at least five annotated ORFs" in result.stderr
    assert not marker.exists()
    assert not output.exists()
    assert not plot.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_refuses_to_overwrite_its_inputs(tmp_path):
    data = tmp_path / "data.csv"
    engine = tmp_path / "estimator.R"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    write_curve(data)
    engine.write_text(
        "get_cutoff_values <- function(path, dest) stop('must not run')\n"
    )
    originals = {path: path.read_bytes() for path in (data, engine)}

    for protected in (data, engine):
        output_alias = run_wrapper(data, protected, plot, engine)
        assert output_alias.returncode != 0
        assert "must not overwrite an input or estimator" in output_alias.stderr

        plot_alias = run_wrapper(data, output, protected, engine)
        assert plot_alias.returncode != 0
        assert "must not overwrite an input or estimator" in plot_alias.stderr

    assert {path: path.read_bytes() for path in (data, engine)} == originals
    assert not output.exists()
    assert not plot.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_rejects_a_preexisting_output_symlink(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "estimator must not load.R"
    engine_marker = tmp_path / "engine-loaded"
    external = tmp_path / "external-parameters.txt"
    write_curve(data)
    external.write_text("0.125,0.4\n")
    output.symlink_to(external)
    plot.write_bytes(VALID_PNG)
    engine.write_text(
        f'writeLines("loaded", {str(engine_marker)!r})\n'
        "get_cutoff_values <- function(path, dest) stop('must not run')\n"
    )

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "output path is not a regular file" in result.stderr
    assert output.is_symlink()
    assert external.read_text() == "0.125,0.4\n"
    assert plot.read_bytes() == VALID_PNG
    assert not engine_marker.exists()
    assert not transaction_path(output, "staging").exists()
    assert not transaction_path(output, "committed").exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_failure_preserves_previous_pair(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "failing estimator.R"
    write_curve(data)
    output.write_bytes(b"old parameters\n")
    plot.write_bytes(b"old plot\n")
    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  writeLines("partial plot", paste0(dest, ".png"))
  stop("simulated nonlinear-fit failure")
}
"""
    )

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "simulated nonlinear-fit failure" in result.stderr
    assert output.read_bytes() == b"old parameters\n"
    assert plot.read_bytes() == b"old plot\n"
    assert list(tmp_path.glob(".parameters.txt.*")) == []
    assert list(tmp_path.glob(".s_curve.png.*")) == []

    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  writeLines("this is definitely not a PNG artifact", paste0(dest, ".png"))
  list(min_RPKM=0.25, min_coverage=0.5)
}
"""
    )
    invalid_plot = run_wrapper(data, output, plot, engine)
    assert invalid_plot.returncode != 0
    assert "is not a PNG" in invalid_plot.stderr
    assert output.read_bytes() == b"old parameters\n"
    assert plot.read_bytes() == b"old plot\n"

    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  writeBin(as.raw(c(
    137,80,78,71,13,10,26,10,
    0,0,0,0,73,69,78,68,174,66,96,130
  )), paste0(dest, ".png"))
  list(min_RPKM=0.25, min_coverage=0.5)
}
"""
    )
    chunkless_plot = run_wrapper(data, output, plot, engine)
    assert chunkless_plot.returncode != 0
    assert "is not a PNG" in chunkless_plot.stderr
    assert output.read_bytes() == b"old parameters\n"
    assert plot.read_bytes() == b"old plot\n"

    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  writeLines("not published", paste0(dest, ".png"))
  list(min_RPKM=Inf, min_coverage=0.5)
}
"""
    )
    invalid_cutoff = run_wrapper(data, output, plot, engine)
    assert invalid_cutoff.returncode != 0
    assert "returned an invalid min_RPKM" in invalid_cutoff.stderr
    assert output.read_bytes() == b"old parameters\n"
    assert plot.read_bytes() == b"old plot\n"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_rejects_a_png_with_a_corrupt_chunk_checksum(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "corrupt PNG estimator.R"
    corrupt_png = tmp_path / "corrupt.png"
    write_curve(data)
    output.write_bytes(b"old parameters\n")
    plot.write_bytes(b"old plot\n")
    corrupted = bytearray(VALID_PNG)
    payload_start = corrupted.index(b"IDAT") + len(b"IDAT")
    corrupted[payload_start] ^= 1
    corrupt_png.write_bytes(corrupted)
    engine.write_text(
        f"""get_cutoff_values <- function(path, dest) {{
  file.copy({str(corrupt_png)!r}, paste0(dest, ".png"))
  list(min_RPKM=0.25, min_coverage=0.5)
}}
"""
    )

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "invalid PNG checksum" in result.stderr
    assert output.read_bytes() == b"old parameters\n"
    assert plot.read_bytes() == b"old plot\n"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
@pytest.mark.parametrize(
    "interruption",
    [
        "first_backup",
        "both_backups",
        "first_install",
        "both_installs",
    ],
)
def test_parameter_wrapper_recovers_interrupted_staging_transaction(
    tmp_path, interruption
):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "failing estimator.R"
    write_curve(data)
    write_failing_engine(engine)

    old_parameters = b"0.125,0.4\n"
    old_plot = VALID_PNG
    output.write_bytes(old_parameters)
    plot.write_bytes(old_plot)
    staging = transaction_path(output, "staging")
    initialize_transaction(staging, (output, plot))

    output.replace(transaction_backup(output, 1))
    if interruption != "first_backup":
        plot.replace(transaction_backup(plot, 2))
    if interruption in {"first_install", "both_installs"}:
        output.write_text("0.25,0.5\n")
    if interruption == "both_installs":
        plot.write_bytes(VALID_PNG)

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "failure after startup reconciliation" in result.stderr
    assert output.read_bytes() == old_parameters
    assert plot.read_bytes() == old_plot
    assert not staging.exists()
    assert not transaction_path(output, "committed").exists()
    assert not transaction_backup(output, 1).exists()
    assert not transaction_backup(plot, 2).exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
@pytest.mark.parametrize("missing", ["data", "engine"])
def test_parameter_wrapper_recovers_before_checking_current_dependencies(
    tmp_path, missing
):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "estimator.R"
    write_curve(data)
    write_failing_engine(engine)
    output.write_text("0.125,0.4\n")
    plot.write_bytes(VALID_PNG)
    staging = transaction_path(output, "staging")
    initialize_transaction(staging, (output, plot))
    output.replace(transaction_backup(output, 1))
    plot.replace(transaction_backup(plot, 2))
    {"data": data, "engine": engine}[missing].unlink()

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "does not exist" in result.stderr
    assert output.read_text() == "0.125,0.4\n"
    assert plot.read_bytes() == VALID_PNG
    assert not staging.exists()
    assert not transaction_backup(output, 1).exists()
    assert not transaction_backup(plot, 2).exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_refuses_a_symlinked_transaction_backup(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "estimator must not run.R"
    marker = tmp_path / "engine-loaded"
    write_curve(data)
    plot.write_bytes(VALID_PNG)
    engine.write_text(
        f'writeLines("loaded", {str(marker)!r})\n'
        "get_cutoff_values <- function(path, dest) stop('must not run')\n"
    )
    external = tmp_path / "external-parameters.txt"
    external.write_text("0.125,0.4\n")
    transaction_backup(output, 1).symlink_to(external)
    staging = transaction_path(output, "staging")
    initialize_transaction(staging, (output, plot))

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "backup is not a regular file" in result.stderr
    assert not marker.exists()
    assert transaction_backup(output, 1).is_symlink()
    assert external.read_text() == "0.125,0.4\n"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_binds_recovery_to_the_original_output_pair(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    original_plot = tmp_path / "original-s-curve.png"
    changed_plot = tmp_path / "changed-s-curve.png"
    engine = tmp_path / "estimator must not run.R"
    marker = tmp_path / "engine-loaded"
    write_curve(data)
    output.write_text("0.125,0.4\n")
    original_plot.write_bytes(VALID_PNG)
    engine.write_text(
        f'writeLines("loaded", {str(marker)!r})\n'
        "get_cutoff_values <- function(path, dest) stop('must not run')\n"
    )
    staging = transaction_path(output, "staging")
    initialize_transaction(staging, (output, original_plot))

    result = run_wrapper(data, output, changed_plot, engine)

    assert result.returncode != 0
    assert "transaction target mismatch for output 2" in result.stderr
    assert output.read_text() == "0.125,0.4\n"
    assert original_plot.read_bytes() == VALID_PNG
    assert not changed_plot.exists()
    assert not marker.exists()
    assert staging.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_refuses_identity_free_retargeted_transaction(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    original_plot = tmp_path / "original-s-curve.png"
    changed_plot = tmp_path / "changed-s-curve.png"
    engine = tmp_path / "estimator must not run.R"
    engine_marker = tmp_path / "engine-loaded"
    write_curve(data)
    output.write_text("0.125,0.4\n")
    original_plot.write_bytes(VALID_PNG)
    changed_plot.write_bytes(CURRENT_PNG)
    engine.write_text(
        f'writeLines("loaded", {str(engine_marker)!r})\n'
        "get_cutoff_values <- function(path, dest) stop('must not run')\n"
    )
    staging = transaction_path(output, "staging")
    staging.mkdir()
    (staging / "present-1").touch()
    (staging / "present-2").touch()
    (staging / "ready").touch()
    output.replace(transaction_backup(output, 1))
    original_plot.replace(transaction_backup(original_plot, 2))

    result = run_wrapper(data, output, changed_plot, engine)

    assert result.returncode != 0
    assert "incomplete target identity" in result.stderr
    assert not engine_marker.exists()
    assert not output.exists()
    assert changed_plot.read_bytes() == CURRENT_PNG
    assert transaction_backup(output, 1).read_text() == "0.125,0.4\n"
    assert transaction_backup(original_plot, 2).read_bytes() == VALID_PNG
    assert staging.is_dir()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_recovers_legacy_random_pair_backups(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "failing estimator.R"
    write_curve(data)
    write_failing_engine(engine)
    parameter_backup = tmp_path / ".parameters.txt.previous.legacy"
    plot_backup = tmp_path / ".s_curve.png.previous.legacy"
    parameter_backup.write_text("0.125,0.4\n")
    plot_backup.write_bytes(VALID_PNG)
    output.write_text("0.25,0.5\n")

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "failure after startup reconciliation" in result.stderr
    assert output.read_text() == "0.125,0.4\n"
    assert plot.read_bytes() == VALID_PNG
    assert not parameter_backup.exists()
    assert not plot_backup.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_can_retry_between_legacy_pair_restores(tmp_path):
    parameter_dir = tmp_path / "parameters"
    plot_dir = tmp_path / "plots"
    parameter_dir.mkdir()
    plot_dir.mkdir()
    data = tmp_path / "data.csv"
    output = parameter_dir / "parameters.txt"
    plot = plot_dir / "s_curve.png"
    engine = tmp_path / "failing estimator.R"
    write_curve(data)
    write_failing_engine(engine)
    parameter_backup = parameter_dir / ".parameters.txt.previous.legacy"
    plot_backup = plot_dir / ".s_curve.png.previous.legacy"
    parameter_backup.write_text("0.125,0.4\n")
    plot_backup.write_bytes(VALID_PNG)
    plot.write_bytes(CURRENT_PNG)

    # Restoration creates the previously missing first target, then cannot
    # replace the second. Both current files are now individually valid, so a
    # durable recovery marker is required to distinguish this mixed pair from a
    # successfully committed publication after a process/power loss.
    plot_dir.chmod(0o500)
    try:
        interrupted = run_wrapper(data, output, plot, engine)
    finally:
        plot_dir.chmod(0o700)

    assert interrupted.returncode != 0
    assert "cannot remove partially published legacy output" in interrupted.stderr
    assert output.read_text() == "0.125,0.4\n"
    assert plot.read_bytes() == CURRENT_PNG
    assert parameter_backup.read_text() == "0.125,0.4\n"
    assert plot_backup.read_bytes() == VALID_PNG
    recovery = parameter_dir / ".parameters.txt.pair-legacy-restore"
    assert recovery.is_file()

    retry = run_wrapper(data, output, plot, engine)

    assert retry.returncode != 0
    assert "failure after startup reconciliation" in retry.stderr
    assert output.read_text() == "0.125,0.4\n"
    assert plot.read_bytes() == VALID_PNG
    assert not recovery.exists()
    assert not parameter_backup.exists()
    assert not plot_backup.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_retries_relative_legacy_pair_restoration(tmp_path):
    workdir = tmp_path / "work"
    parameter_dir = workdir / "parameters"
    plot_dir = workdir / "plots"
    parameter_dir.mkdir(parents=True)
    plot_dir.mkdir()
    data = workdir / "data.csv"
    engine = workdir / "failing estimator.R"
    write_curve(data)
    write_failing_engine(engine)
    output = Path("parameters/parameters.txt")
    plot = Path("plots/s_curve.png")
    parameter_backup = parameter_dir / ".parameters.txt.previous.legacy"
    plot_backup = plot_dir / ".s_curve.png.previous.legacy"
    parameter_backup.write_text("0.125,0.4\n")
    plot_backup.write_bytes(VALID_PNG)
    (workdir / plot).write_bytes(CURRENT_PNG)

    plot_dir.chmod(0o500)
    try:
        interrupted = run_wrapper(
            Path("data.csv"), output, plot, Path("failing estimator.R"), cwd=workdir
        )
    finally:
        plot_dir.chmod(0o700)

    assert interrupted.returncode != 0
    assert "cannot remove partially published legacy output" in interrupted.stderr
    recovery = parameter_dir / ".parameters.txt.pair-legacy-restore"
    assert recovery.read_text().splitlines() == [
        str((workdir / output).resolve()),
        str((workdir / plot).resolve()),
    ]
    assert parameter_backup.exists()
    assert plot_backup.exists()

    retry = run_wrapper(
        Path("data.csv"), output, plot, Path("failing estimator.R"), cwd=workdir
    )

    assert retry.returncode != 0
    assert "failure after startup reconciliation" in retry.stderr
    assert (workdir / output).read_text() == "0.125,0.4\n"
    assert (workdir / plot).read_bytes() == VALID_PNG
    assert not recovery.exists()
    assert not parameter_backup.exists()
    assert not plot_backup.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_refuses_an_incomplete_legacy_snapshot(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "estimator must not load.R"
    marker = tmp_path / "engine-loaded"
    legacy_output = tmp_path / ".parameters.txt.previous.interrupted"
    write_curve(data)
    plot.write_bytes(VALID_PNG)
    legacy_output.write_text("0.125,0.4\n")
    engine.write_text(
        f'writeLines("loaded", {str(marker)!r})\n'
        "get_cutoff_values <- function(path, dest) stop('must not run')\n"
    )

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "legacy publication snapshot is incomplete" in result.stderr
    assert not output.exists()
    assert plot.read_bytes() == VALID_PNG
    assert legacy_output.read_text() == "0.125,0.4\n"
    assert not marker.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_fails_closed_on_an_orphaned_deterministic_backup(
    tmp_path,
):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "estimator must not run.R"
    marker = tmp_path / "engine-loaded"
    write_curve(data)
    plot.write_bytes(VALID_PNG)
    orphan = transaction_backup(output, 1)
    orphan.write_text("0.125,0.4\n")
    engine.write_text(
        f'writeLines("loaded", {str(marker)!r})\n'
        "get_cutoff_values <- function(path, dest) stop('must not run')\n"
    )

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "backup exists without transaction state" in result.stderr
    assert orphan.read_text() == "0.125,0.4\n"
    assert plot.read_bytes() == VALID_PNG
    assert not output.exists()
    assert not marker.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_serializes_concurrent_publishers(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "slow estimator.R"
    activity = tmp_path / "estimator-activity.txt"
    write_curve(data)
    engine.write_text(
        f"""get_cutoff_values <- function(path, dest) {{
  cat("start\\n", file={str(activity)!r}, append=TRUE)
  Sys.sleep(0.75)
  cat("end\\n", file={str(activity)!r}, append=TRUE)
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.125, min_coverage=0.4)
}}
"""
    )
    command = wrapper_command(data, output, plot, engine)

    first = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    second = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    first_stdout, first_stderr = first.communicate(timeout=20)
    second_stdout, second_stderr = second.communicate(timeout=20)

    assert first.returncode == 0, (first_stdout, first_stderr)
    assert second.returncode == 0, (second_stdout, second_stderr)
    assert activity.read_text().splitlines() == ["start", "end", "start", "end"]
    assert output.read_text() == "0.125,0.40000000000000002\n"
    assert plot.read_bytes().startswith(b"\x89PNG")
    assert not transaction_path(output, "staging").exists()
    assert not transaction_path(output, "committed").exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_launcher_locks_tilde_and_absolute_output_aliases_together(
    tmp_path,
):
    home = tmp_path / "home"
    workdir = tmp_path / "work"
    home.mkdir()
    workdir.mkdir()
    data = tmp_path / "data.csv"
    engine = tmp_path / "slow estimator.R"
    activity = tmp_path / "estimator-activity.txt"
    write_curve(data)
    engine.write_text(
        f"""get_cutoff_values <- function(path, dest) {{
  cat("start\\n", file={str(activity)!r}, append=TRUE)
  Sys.sleep(0.75)
  cat("end\\n", file={str(activity)!r}, append=TRUE)
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.125, min_coverage=0.4)
}}
"""
    )
    absolute_directory = home / "cutoffs"
    environment = {**os.environ, "HOME": str(home)}
    environment.pop("R_USER", None)
    tilde_command = wrapper_command(
        data,
        "~/cutoffs/parameters.txt",
        "~/cutoffs/s_curve.png",
        engine,
    )
    absolute_command = wrapper_command(
        data,
        absolute_directory / "parameters.txt",
        absolute_directory / "s_curve.png",
        engine,
    )

    first = subprocess.Popen(
        tilde_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=workdir,
        env=environment,
    )
    deadline = time.monotonic() + 10
    while not activity.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert activity.exists(), first.communicate(timeout=1)
    second = subprocess.Popen(
        absolute_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=workdir,
        env=environment,
    )
    time.sleep(0.2)
    assert activity.read_text().splitlines() == ["start"]

    first_stdout, first_stderr = first.communicate(timeout=10)
    second_stdout, second_stderr = second.communicate(timeout=10)

    assert first.returncode == 0, (first_stdout, first_stderr)
    assert second.returncode == 0, (second_stdout, second_stderr)
    assert activity.read_text().splitlines() == ["start", "end", "start", "end"]
    assert (absolute_directory / "parameters.txt").is_file()
    assert (absolute_directory / "s_curve.png").is_file()
    assert not (workdir / "~").exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_launcher_releases_its_lock_when_the_estimator_is_killed(
    tmp_path,
):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    slow_engine = tmp_path / "slow estimator.R"
    fast_engine = tmp_path / "fast estimator.R"
    started = tmp_path / "started"
    write_curve(data)
    slow_engine.write_text(
        f"""get_cutoff_values <- function(path, dest) {{
  writeLines("started", {str(started)!r})
  Sys.sleep(30)
  stop("the killed estimator unexpectedly resumed")
}}
"""
    )
    fast_engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.25, min_coverage=0.5)
}
"""
    )

    first = subprocess.Popen(
        wrapper_command(data, output, plot, slow_engine),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    deadline = time.monotonic() + 10
    while not started.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert started.exists(), first.communicate(timeout=1)

    second = subprocess.Popen(
        wrapper_command(data, output, plot, fast_engine),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    time.sleep(0.2)
    assert second.poll() is None

    first.kill()
    first.communicate(timeout=5)
    second_stdout, second_stderr = second.communicate(timeout=10)

    assert first.returncode < 0
    assert second.returncode == 0, (second_stdout, second_stderr)
    assert output.read_text() == "0.25,0.5\n"
    assert plot.read_bytes().startswith(b"\x89PNG")


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_removes_partial_outputs_that_were_new(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "failing estimator.R"
    write_curve(data)
    write_failing_engine(engine)
    staging = transaction_path(output, "staging")
    initialize_transaction(staging, (output, plot), existed=(False, False))
    output.write_text("0.25,0.5\n")

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "failure after startup reconciliation" in result.stderr
    assert not output.exists()
    assert not plot.exists()
    assert not staging.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
@pytest.mark.parametrize("remaining_backups", [(1, 2), (2,)])
def test_parameter_wrapper_keeps_a_committed_pair_and_cleans_backups(
    tmp_path, remaining_backups
):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "failing estimator.R"
    write_curve(data)
    write_failing_engine(engine)

    output.write_text("0.25,0.5\n")
    plot.write_bytes(VALID_PNG)
    committed = transaction_path(output, "committed")
    initialize_transaction(committed, (output, plot))
    if 1 in remaining_backups:
        transaction_backup(output, 1).write_text("0.125,0.4\n")
    if 2 in remaining_backups:
        transaction_backup(plot, 2).write_bytes(VALID_PNG)

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "failure after startup reconciliation" in result.stderr
    assert output.read_text() == "0.25,0.5\n"
    assert plot.read_bytes() == VALID_PNG
    assert not committed.exists()
    assert not transaction_backup(output, 1).exists()
    assert not transaction_backup(plot, 2).exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_defers_cleanup_after_a_valid_commit(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "cleanup-fault estimator.R"
    write_curve(data)
    output.write_text("0.125,0.4\n")
    plot.write_bytes(VALID_PNG)
    engine.write_text(
        """base_unlink <- base::unlink
assign("unlink", function(x, recursive=FALSE, force=FALSE) {
  if (any(grepl(".pair-transaction-old-", x, fixed=TRUE))) return(1L)
  base_unlink(x, recursive=recursive, force=force)
}, envir=.GlobalEnv)
get_cutoff_values <- function(path, dest) {
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.25, min_coverage=0.5)
}
"""
    )

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode == 0, result.stderr
    assert "committed outputs are valid" in result.stderr
    assert output.read_text() == "0.25,0.5\n"
    assert plot.read_bytes().startswith(b"\x89PNG")
    assert transaction_path(output, "committed").is_dir()
    assert transaction_backup(output, 1).exists()
    assert transaction_backup(plot, 2).exists()

    engine.write_text(
        """get_cutoff_values <- function(path, dest) {
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.375, min_coverage=0.6)
}
"""
    )
    retry = run_wrapper(data, output, plot, engine)

    assert retry.returncode == 0, retry.stderr
    assert output.read_text() == "0.375,0.59999999999999998\n"
    assert not transaction_path(output, "committed").exists()
    assert not transaction_backup(output, 1).exists()
    assert not transaction_backup(plot, 2).exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_recovers_interrupted_transaction_metadata_cleanup(
    tmp_path,
):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "cleanup-fault estimator.R"
    write_curve(data)
    engine.write_text(
        """base_unlink <- base::unlink
assign("unlink", function(x, recursive=FALSE, force=FALSE) {
  if (recursive && any(grepl(".pair-transaction.committed", x, fixed=TRUE))) {
    return(1L)
  }
  base_unlink(x, recursive=recursive, force=force)
}, envir=.GlobalEnv)
get_cutoff_values <- function(path, dest) {
  png(paste0(dest, ".png"), width=600, height=600)
  plot(1, 1)
  invisible(dev.off())
  list(min_RPKM=0.25, min_coverage=0.5)
}
"""
    )

    interrupted = run_wrapper(data, output, plot, engine)

    assert interrupted.returncode == 0, interrupted.stderr
    assert "transaction cleanup is deferred" in interrupted.stderr
    committed = transaction_path(output, "committed")
    assert committed.is_dir()
    assert not (committed / "ready").exists()

    # Model SIGKILL partway through recursive deletion: the durable no-ready
    # state can have any subset of its now-disposable marker files left.
    (committed / "absent-1").unlink()
    (committed / "target-2").unlink()
    write_failing_engine(engine)

    retry = run_wrapper(data, output, plot, engine)

    assert retry.returncode != 0
    assert "failure after startup reconciliation" in retry.stderr
    assert "incomplete state" not in retry.stderr
    assert output.read_text() == "0.25,0.5\n"
    assert plot.read_bytes().startswith(b"\x89PNG")
    assert not committed.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_restores_snapshot_if_committed_pair_is_invalid(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "failing estimator.R"
    write_curve(data)
    write_failing_engine(engine)

    old_parameters = b"0.125,0.4\n"
    old_plot = VALID_PNG
    output.write_text("truncated parameter artifact\n")
    plot.write_bytes(VALID_PNG)
    committed = transaction_path(output, "committed")
    initialize_transaction(committed, (output, plot))
    transaction_backup(output, 1).write_bytes(old_parameters)
    transaction_backup(plot, 2).write_bytes(old_plot)

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "failure after startup reconciliation" in result.stderr
    assert output.read_bytes() == old_parameters
    assert plot.read_bytes() == old_plot
    assert not committed.exists()
    assert not transaction_backup(output, 1).exists()
    assert not transaction_backup(plot, 2).exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_discards_interrupted_transaction_setup(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "failing estimator.R"
    write_curve(data)
    write_failing_engine(engine)
    output.write_text("0.125,0.4\n")
    plot.write_bytes(VALID_PNG)
    staging = transaction_path(output, "staging")
    staging.mkdir()
    (staging / "present-1").touch()

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "failure after startup reconciliation" in result.stderr
    assert output.read_text() == "0.125,0.4\n"
    assert plot.read_bytes() == VALID_PNG
    assert not staging.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_parameter_wrapper_fails_closed_on_ambiguous_transaction(tmp_path):
    data = tmp_path / "data.csv"
    output = tmp_path / "parameters.txt"
    plot = tmp_path / "s_curve.png"
    engine = tmp_path / "estimator must not run.R"
    marker = tmp_path / "engine-loaded"
    write_curve(data)
    output.write_text("0.125,0.4\n")
    plot.write_bytes(VALID_PNG)
    engine.write_text(
        f'writeLines("loaded", {str(marker)!r})\n'
        "get_cutoff_values <- function(path, dest) stop('must not run')\n"
    )
    staging = transaction_path(output, "staging")
    initialize_transaction(staging, (output, plot))
    output.unlink()

    result = run_wrapper(data, output, plot, engine)

    assert result.returncode != 0
    assert "publication transaction lost existing output 1" in result.stderr
    assert not marker.exists()
    assert not output.exists()
    assert plot.read_bytes() == VALID_PNG
    assert staging.exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is required")
def test_production_parameter_rule_recovers_after_postcommit_process_group_kill(
    tmp_path, snakemake_command
):
    invocation, run_environment, workflow_root, receipt, snakefile = (
        production_parameter_rule_run(tmp_path, snakemake_command)
    )
    successful = subprocess.run(
        invocation,
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_PARAMETER_MODE": "first"},
        timeout=120,
    )
    assert successful.returncode == 0, successful.stdout + successful.stderr

    cutoff_directory = receipt.parent
    parameters = cutoff_directory / "parameters.txt"
    plot = cutoff_directory / "s_curve.png"
    old_only = cutoff_directory / "old-only.txt"
    old_only.write_text("preserve runner-managed siblings\n")
    old_snapshot = {
        "parameters": parameters.read_bytes(),
        "plot": plot.read_bytes(),
        "old_only": old_only.read_bytes(),
    }

    notification = workflow_root / "parameter-runner-committed.pid"
    process = subprocess.Popen(
        [*invocation, "--forcerun", "parameterEstimation"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=workflow_root,
        env={
            **run_environment,
            "FAKE_PARAMETER_MODE": "second",
            "HRIBO_TEST_PARAMETER_INSTALL_NOTIFY": str(notification),
        },
        start_new_session=True,
    )
    try:
        deadline = time.monotonic() + 60
        while not notification.exists() and process.poll() is None:
            if time.monotonic() >= deadline:
                break
            time.sleep(0.01)
        assert notification.exists(), "runner never published its completion receipt"
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate(timeout=30)
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait(timeout=10)

    assert process.returncode < 0, stdout + stderr
    assert parameters.read_text() == "0.25,0.5\n"
    assert plot.read_bytes().startswith(b"\x89PNG")
    assert plot.read_bytes() != old_snapshot["plot"]
    assert old_only.read_bytes() == old_snapshot["old_only"]
    assert not transaction_path(parameters, "staging").exists()
    assert not transaction_path(parameters, "committed").exists()

    unlocked = subprocess.run(
        [*snakemake_command, "--snakefile", str(snakefile), "--unlock"],
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env=run_environment,
        timeout=30,
    )
    assert unlocked.returncode == 0, unlocked.stdout + unlocked.stderr

    retry = subprocess.run(
        [*invocation, "--forcerun", "parameterEstimation"],
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_PARAMETER_MODE": "fail"},
        timeout=120,
    )

    assert retry.returncode != 0
    assert "simulated parameter-estimation failure" in retry.stderr
    assert parameters.read_text() == "0.25,0.5\n"
    assert plot.read_bytes().startswith(b"\x89PNG")
    assert old_only.read_bytes() == old_snapshot["old_only"]
    assert not receipt.exists()
    backup_root = workflow_root / ".snakemake/backups"
    assert not backup_root.exists() or not any(backup_root.rglob("*"))


def test_parameter_rule_tracks_patched_engine_and_diagnostic_plot():
    assert "rule prepareDeepRiboSCurveScript:" in RULE
    assert (
        'patcher=workflow.source_path("../scripts/patch_deepribo_scurve.py")'
        in RULE
    )
    assert "/usr/local/bin/s_curve_cutoff_estimation.R" in RULE
    assert "engine=rules.prepareDeepRiboSCurveScript.output.script" in RULE
    assert (
        'launcher=workflow.source_path("../scripts/run_parameter_estimation.py")'
        in RULE
    )
    assert (
        'wrapper_deps=[workflow.source_path("../scripts/parameter_estimation.R")]'
        in RULE
    )
    assert "python3 {input.launcher:q} {input.wrapper_deps:q}" in RULE
    assert 'receipt=ensure(' in RULE
    assert '"deepribo/cutoffs/{condition}-{replicate}/.complete"' in RULE
    assert "non_empty=True" in RULE
    assert 'os.path.dirname(output.receipt), "parameters.txt"' in RULE
    assert 'os.path.dirname(output.receipt), "s_curve.png"' in RULE
    assert "--out {params.parameters:q}" in RULE
    assert "--plot {params.plot:q}" in RULE
    assert "--receipt {output.receipt:q}" in RULE
    assert "{input.engine:q}" in RULE
    assert "cutoff_receipt=rules.parameterEstimation.output.receipt" in RULE
    assert "input.cutoff_receipt" in RULE
    assert 'update("deepribo/' not in RULE


def test_parameter_transaction_persists_state_before_ready_barrier():
    script = WRAPPER.read_text()
    publication = script[script.index("publish_pair <- function"):]
    state_fsync = publication.index("files = transaction_state_files")
    ready_create = publication.index("create_transaction_marker(ready)")
    ready_fsync = publication.index("fsync_paths(files = ready")
    first_output_mutation = publication.index(
        "transaction_backup(target, index)", ready_fsync
    )

    assert state_fsync < ready_create < ready_fsync < first_output_mutation
