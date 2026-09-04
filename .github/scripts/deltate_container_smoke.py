#!/usr/bin/env python3
"""Exercise HRIBO's production DeltaTE rules through Apptainer."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"
CONTRAST = "A-B"
NEUTRAL_IDS = tuple(f"neutral_{index:03d}" for index in range(48))
UP_IDS = tuple(f"up_{index:03d}" for index in range(16))
DOWN_IDS = tuple(f"down_{index:03d}" for index in range(16))
MATCHED_UP_IDS = tuple(f"matched_up_{index:03d}" for index in range(8))
MATCHED_DOWN_IDS = tuple(f"matched_down_{index:03d}" for index in range(8))
MATCHED_IDS = (*MATCHED_UP_IDS, *MATCHED_DOWN_IDS)
LOW_COUNT_IDS = ("low_count_000",)
EFFECT_IDS = (*NEUTRAL_IDS, *UP_IDS, *DOWN_IDS, *MATCHED_IDS)
ALL_IDS = {*EFFECT_IDS, *LOW_COUNT_IDS}
EXPECTED_RESULT_ROWS = 97


def _negative_binomial_count(rng: np.random.Generator, mean: float) -> int:
    """Draw one stable, overdispersed count for the bounded fixture."""

    shape = 80.0
    probability = shape / (shape + mean)
    return max(2, int(rng.negative_binomial(shape, probability)))


def _write_count_fixture(path: Path) -> None:
    """Write balanced positive, negative, and neutral translation effects."""

    rng = np.random.default_rng(20260903)
    rows: list[dict[str, int | str]] = []
    replicate_scales = {"1": 0.97, "2": 1.03}

    for gene_index, identifier in enumerate(EFFECT_IDS):
        if identifier.startswith("neutral_"):
            paired_index = gene_index
            a_ribo_effect = b_ribo_effect = 1.0
            a_rna_effect = b_rna_effect = 1.0
        elif identifier.startswith("up_"):
            paired_index = gene_index - len(NEUTRAL_IDS)
            a_ribo_effect, b_ribo_effect = 4.0, 1.0
            a_rna_effect = b_rna_effect = 1.0
        elif identifier.startswith("down_"):
            paired_index = gene_index - len(NEUTRAL_IDS) - len(UP_IDS)
            a_ribo_effect, b_ribo_effect = 1.0, 4.0
            a_rna_effect = b_rna_effect = 1.0
        elif identifier.startswith("matched_up_"):
            paired_index = (
                gene_index - len(NEUTRAL_IDS) - len(UP_IDS) - len(DOWN_IDS)
            )
            a_ribo_effect, b_ribo_effect = 4.0, 1.0
            a_rna_effect, b_rna_effect = 4.0, 1.0
        else:
            paired_index = (
                gene_index
                - len(NEUTRAL_IDS)
                - len(UP_IDS)
                - len(DOWN_IDS)
                - len(MATCHED_UP_IDS)
            )
            a_ribo_effect, b_ribo_effect = 1.0, 4.0
            a_rna_effect, b_rna_effect = 1.0, 4.0

        rna_mean = 180.0 + 11.0 * (paired_index % 17)
        ribo_mean = 150.0 + 13.0 * (paired_index % 19)
        row: dict[str, int | str] = {"Identifier": identifier}
        condition_effects = (
            ("A", a_ribo_effect, a_rna_effect),
            ("B", b_ribo_effect, b_rna_effect),
        )
        for condition, ribo_effect, rna_effect in condition_effects:
            for replicate, scale in replicate_scales.items():
                row[f"RIBO-{condition}-{replicate}"] = _negative_binomial_count(
                    rng, ribo_mean * ribo_effect * scale
                )
                row[f"RNA-{condition}-{replicate}"] = _negative_binomial_count(
                    rng, rna_mean * rna_effect * scale
                )
        rows.append(row)

    rows.append(
        {
            "Identifier": LOW_COUNT_IDS[0],
            "RIBO-A-1": 4,
            "RNA-A-1": 5,
            "RIBO-A-2": 6,
            "RNA-A-2": 7,
            "RIBO-B-1": 4,
            "RNA-B-1": 5,
            "RIBO-B-2": 6,
            "RNA-B-2": 7,
        }
    )

    if len(rows) != EXPECTED_RESULT_ROWS:
        raise AssertionError(
            f"fixture has {len(rows)} rows, expected {EXPECTED_RESULT_ROWS}"
        )
    path.parent.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def _write_inputs(workdir: Path) -> Path:
    """Create all inputs without invoking mapping or external downloads."""

    inputs = workdir / "input data"
    fastq_dir = inputs / "fastq"
    bam_dir = workdir / "bam"
    fastq_dir.mkdir(parents=True)
    bam_dir.mkdir(parents=True)

    genome = inputs / "tiny genome.fa"
    genome.write_text(
        ">fixture_contig deltaTE smoke fixture\n" + "ATG" * 200 + "\n"
    )
    annotation = inputs / "tiny annotation.gff"
    annotation.write_text(
        "##gff-version 3\n"
        "fixture_contig\tfixture\tCDS\t1\t300\t.\t+\t0\t"
        "ID=fixture_cds;locus_tag=fixture_cds\n"
    )

    sample_rows: list[dict[str, str]] = []
    for method in ("RIBO", "RNA"):
        for condition in ("A", "B"):
            for replicate in ("1", "2"):
                stem = f"{method}-{condition}-{replicate}"
                fastq = fastq_dir / f"{stem}.fastq.gz"
                read = f"@{stem}\n{'A' * 30}\n+\n{'I' * 30}\n".encode()
                fastq.write_bytes(gzip.compress(read, mtime=0))
                sample_rows.append(
                    {
                        "method": method,
                        "condition": condition,
                        "replicate": replicate,
                        "fastqFile": str(fastq),
                        "fastqFile2": "",
                    }
                )
                # The adapter uses declared BAM basenames, not BAM contents.
                (bam_dir / f"{stem}.bam").touch()

    samples = inputs / "samples.tsv"
    pd.DataFrame(sample_rows).to_csv(samples, sep="\t", index=False)

    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["biologySettings"].update(
        {
            "genome": str(genome),
            "annotation": str(annotation),
            "samples": str(samples),
        }
    )
    config["differentialExpressionSettings"].update(
        {"contrasts": [CONTRAST], "xtailBins": 1000}
    )
    config["predictionSettings"]["deepribo"] = "off"
    config["workflowSettings"]["stages"] = ["differential_expression"]
    config_path = inputs / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    _write_count_fixture(
        workdir / "readcounts" / "differential_expression_read_counts.csv"
    )
    return config_path


def _stream(command: list[str], environment: dict[str, str]) -> str:
    """Run a command while retaining its combined output for assertions."""

    print("+ " + " ".join(command), flush=True)
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=environment,
    )
    assert process.stdout is not None
    rendered: list[str] = []
    for line in process.stdout:
        print(line, end="", flush=True)
        rendered.append(line)
    returncode = process.wait()
    output = "".join(rendered)
    if returncode != 0:
        raise RuntimeError(f"command failed with exit status {returncode}")
    return output


def _read_result(path: Path, header: list[str]) -> dict[str, list[float]]:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != header:
        raise AssertionError(f"unexpected header in {path}: {rows[:1]!r}")
    if len(rows) != EXPECTED_RESULT_ROWS + 1:
        raise AssertionError(
            f"{path} has {len(rows) - 1} rows, expected {EXPECTED_RESULT_ROWS}"
        )

    values: dict[str, list[float]] = {}
    for row in rows[1:]:
        if len(row) != len(header) + 1:
            raise AssertionError(f"malformed result row in {path}: {row!r}")
        identifier = row[0]
        numbers: list[float] = []
        for column, value in zip(header, row[1:], strict=True):
            # DESeq2 legitimately leaves adjusted p-values undefined for rows
            # removed by independent filtering. All fitted statistics and raw
            # p-values must still be present and finite.
            if column == "padj" and value == "NA":
                numbers.append(math.nan)
                continue
            number = float(value)
            if not math.isfinite(number):
                raise AssertionError(
                    f"non-finite {column} for {identifier} in {path}"
                )
            if column in {"pvalue", "padj"} and not 0 <= number <= 1:
                raise AssertionError(
                    f"out-of-range {column} for {identifier} in {path}"
                )
            numbers.append(number)
        if identifier in values:
            raise AssertionError(f"duplicate result identifier {identifier!r}")
        values[identifier] = numbers
    return values


def _outputs(workdir: Path) -> tuple[Path, ...]:
    return (
        workdir / "deltate/A-B/fold_changes/deltaRibo.txt",
        workdir / "deltate/A-B/fold_changes/deltaRNA.txt",
        workdir / "deltate/A-B/fold_changes/deltaTE.txt",
        workdir / "deltate/A-B_figures.pdf",
    )


def _validate_outputs(workdir: Path) -> None:
    outputs = _outputs(workdir)
    common_header = ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]
    ribo = _read_result(outputs[0], common_header)
    rna = _read_result(outputs[1], common_header)
    te = _read_result(
        outputs[2],
        ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"],
    )
    for name, table in (("deltaRibo", ribo), ("deltaRNA", rna), ("deltaTE", te)):
        if set(table) != ALL_IDS:
            raise AssertionError(f"{name} identifiers differ from the fixture")

    up_median = statistics.median(te[name][1] for name in UP_IDS)
    down_median = statistics.median(te[name][1] for name in DOWN_IDS)
    neutral_median = statistics.median(te[name][1] for name in NEUTRAL_IDS)
    if up_median <= 0:
        raise AssertionError(f"A-up deltaTE effects are not positive: {up_median}")
    if down_median >= 0:
        raise AssertionError(f"A-down deltaTE effects are not negative: {down_median}")
    if abs(neutral_median) >= 0.5:
        raise AssertionError(
            f"neutral deltaTE median is unexpectedly large: {neutral_median}"
        )

    # These controls change abundance in the same direction and by the same
    # factor in RNA and RIBO. Both marginal models must recover the change,
    # while their interaction (translation efficiency) must remain near zero.
    for label, identifiers, direction in (
        ("matched-up", MATCHED_UP_IDS, 1),
        ("matched-down", MATCHED_DOWN_IDS, -1),
    ):
        ribo_median = statistics.median(ribo[name][1] for name in identifiers)
        rna_median = statistics.median(rna[name][1] for name in identifiers)
        te_median = statistics.median(te[name][1] for name in identifiers)
        if direction * ribo_median <= 1 or direction * rna_median <= 1:
            raise AssertionError(
                f"{label} RNA/RIBO effects are not strong and concordant: "
                f"RIBO={ribo_median}, RNA={rna_median}"
            )
        if abs(ribo_median - rna_median) >= 0.5:
            raise AssertionError(
                f"{label} RNA/RIBO effects are not matched: "
                f"RIBO={ribo_median}, RNA={rna_median}"
            )
        if abs(te_median) >= 0.5:
            raise AssertionError(
                f"{label} deltaTE median is unexpectedly large: {te_median}"
            )

    engine = workdir / "deltate" / "DTEG.R"
    engine_text = engine.read_text()
    for marker in (
        "varianceStabilizingTransformation(ddsMat_ribo, blind=TRUE)",
        "varianceStabilizingTransformation(ddsMat_rna, blind=TRUE)",
        'plot_class_example(forwarded, "Forwarded")',
        "max(abs(c(res_ribo[,2],res_rna[,2]))",
    ):
        if marker not in engine_text:
            raise AssertionError(f"patched DeltaTE engine lacks marker: {marker}")

    pdf = outputs[3].read_bytes()
    if not pdf.startswith(b"%PDF-") or b"%%EOF" not in pdf[-2048:]:
        raise AssertionError("published DeltaTE figure is not a complete PDF")
    pdfinfo = shutil.which("pdfinfo")
    if pdfinfo is None:
        raise RuntimeError("pdfinfo is required to validate the DeltaTE report")
    info = subprocess.run(
        [pdfinfo, str(outputs[3])],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    pages = next(
        (
            line.split(":", 1)[1].strip()
            for line in info.splitlines()
            if line.startswith("Pages:")
        ),
        None,
    )
    if pages != "4":
        raise AssertionError(f"DeltaTE figure has {pages!r} pages, expected four")


def _snapshot(workdir: Path) -> dict[str, tuple[str, int]]:
    paths = (*_outputs(workdir), workdir / "deltate" / "DTEG.R")
    return {
        str(path): (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
        for path in paths
    }


def _empty_external_workdir(path: Path) -> Path:
    workdir = path.expanduser().resolve()
    if workdir == REPO or REPO in workdir.parents:
        raise ValueError("the smoke work directory must be outside the checkout")
    if workdir.exists() and any(workdir.iterdir()):
        raise ValueError(f"the smoke work directory is not empty: {workdir}")
    workdir.mkdir(parents=True, exist_ok=True)
    return workdir


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workdir",
        type=Path,
        required=True,
        help="new or empty directory outside the repository checkout",
    )
    parser.add_argument(
        "--conda-prefix",
        type=Path,
        default=REPO / ".snakemake" / "ci-deltate-conda",
    )
    parser.add_argument(
        "--apptainer-prefix",
        type=Path,
        default=REPO / ".snakemake" / "ci-deltate-apptainer",
    )
    parser.add_argument("--snakemake", default="snakemake")
    return parser


def main() -> None:
    args = _parser().parse_args()
    workdir = _empty_external_workdir(args.workdir)
    conda_prefix = args.conda_prefix.expanduser().resolve()
    apptainer_prefix = args.apptainer_prefix.expanduser().resolve()
    snakemake = shutil.which(args.snakemake)
    if snakemake is None:
        raise RuntimeError(f"Snakemake executable is unavailable: {args.snakemake}")
    apptainer = shutil.which("apptainer")
    if apptainer is None:
        raise RuntimeError("Apptainer is unavailable on PATH")

    version = subprocess.run(
        [apptainer, "version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    print(f"Using Apptainer {version}", flush=True)
    config_path = _write_inputs(workdir)
    conda_prefix.mkdir(parents=True, exist_ok=True)
    apptainer_prefix.mkdir(parents=True, exist_ok=True)

    cache_sibling = workdir.parent / f".{workdir.name}-runtime-cache"
    environment = {
        **os.environ,
        "XDG_CACHE_HOME": os.environ.get(
            "XDG_CACHE_HOME", str(cache_sibling / "xdg")
        ),
        "APPTAINER_CACHEDIR": os.environ.get(
            "APPTAINER_CACHEDIR", str(cache_sibling / "apptainer")
        ),
        "APPTAINER_TMPDIR": os.environ.get(
            "APPTAINER_TMPDIR", str(cache_sibling / "apptainer-tmp")
        ),
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    for variable in ("XDG_CACHE_HOME", "APPTAINER_CACHEDIR", "APPTAINER_TMPDIR"):
        Path(environment[variable]).mkdir(parents=True, exist_ok=True)

    targets = [str(path.relative_to(workdir)) for path in _outputs(workdir)]
    command = [
        snakemake,
        *targets,
        "--cores",
        "1",
        "--printshellcmds",
        "--show-failed-logs",
        "--rerun-incomplete",
        "--notemp",
        "--latency-wait",
        "60",
        "--software-deployment-method",
        "conda",
        "apptainer",
        "--conda-prefix",
        str(conda_prefix),
        "--apptainer-prefix",
        str(apptainer_prefix),
        "--allowed-rules",
        "contrastInput",
        "deltatePrepareInput",
        "prepareDeltaTEScript",
        "deltate",
        "--snakefile",
        str(SNAKEFILE),
        "--directory",
        str(workdir),
        "--configfile",
        str(config_path),
    ]

    first = _stream(command, environment)
    for rule in (
        "contrastInput",
        "deltatePrepareInput",
        "prepareDeltaTEScript",
        "deltate",
    ):
        if f"rule {rule}:" not in first:
            raise AssertionError(f"production rule did not execute: {rule}")
    for helper in ("patch_deltate.py", "run_deltate.sh"):
        cached_command = any(
            "source-cache" in line and helper in line for line in first.splitlines()
        )
        if f"{helper} (cached)" not in first or not cached_command:
            raise AssertionError(f"container helper was not source-cached: {helper}")
    _validate_outputs(workdir)
    before = _snapshot(workdir)

    second = _stream(command, environment)
    if "Nothing to be done" not in second:
        raise AssertionError("second Snakemake run was not a no-op")
    if _snapshot(workdir) != before:
        raise AssertionError("a no-op rerun changed a published artifact")

    print("DeltaTE Snakemake + Apptainer smoke passed", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (
        AssertionError,
        OSError,
        RuntimeError,
        ValueError,
        subprocess.SubprocessError,
    ) as error:
        print(f"DeltaTE container smoke failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
