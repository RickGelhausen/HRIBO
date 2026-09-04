#!/usr/bin/env python3
"""Exercise HRIBO's production DeepRibo rules through Apptainer."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import shutil
import subprocess
import sys
import zlib
from dataclasses import dataclass
from pathlib import Path

import pysam
import yaml


REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"
CONTIG = "synthetic_contig"
CONDITION = "A"
REPLICATE = "1"
CODING_LENGTH = 300
SPACER_LENGTH = 60
PROFILE_LENGTHS = (
    13,
    15,
    20,
    28,
    41,
    61,
    89,
    125,
    163,
    199,
    227,
    247,
    260,
    268,
    273,
    275,
)
PROFILE_COUNTS = (3, 4, 4, 5, 5, 5, 5, 5, 6, 7, 10, 13, 19, 27, 40, 59)
PARSED_COLUMNS = (
    "filename",
    "filename_counts",
    "label",
    "in_gene",
    "strand",
    "coverage",
    "coverage_elo",
    "rpk",
    "rpk_elo",
    "start_site",
    "start_codon",
    "stop_site",
    "stop_codon",
    "locus",
    "prot_seq",
    "nuc_seq",
)
PREDICTION_COLUMNS = (
    *PARSED_COLUMNS,
    "pred",
    "pred_rank",
    "SS",
    "dist",
    "SS_pred_rank",
)


@dataclass(frozen=True)
class FixtureOrf:
    start: int
    stop: int
    strand: str
    annotated: bool
    signal_length: int
    signal_count: int

    @property
    def parser_locus(self) -> str:
        if self.strand == "+":
            return f"{CONTIG}:{self.start}-{self.stop - 2}"
        return f"{CONTIG}:{self.start + 2}-{self.stop}"

    @property
    def identifier(self) -> str:
        return f"{CONTIG}:{self.start}-{self.stop}:{self.strand}"

    @property
    def asite_coverage(self) -> float:
        return self.signal_length / CODING_LENGTH

    @property
    def asite_occupancy(self) -> float:
        return self.signal_length * self.signal_count / CODING_LENGTH

    @property
    def read_coverage(self) -> float:
        return (self.signal_length + 12) / CODING_LENGTH

    @property
    def read_occupancy(self) -> float:
        overlap = sum(
            min(position + 13, CODING_LENGTH) - max(position - 17, 0)
            for position in range(self.signal_length)
        )
        return overlap * self.signal_count / CODING_LENGTH


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _build_genome() -> tuple[str, tuple[FixtureOrf, ...]]:
    """Build isolated annotated ORFs with a smooth occupancy/coverage curve."""

    coding_orf = "ATG" + "AAA" * 99 + "TAA"
    pieces = ["N" * SPACER_LENGTH]
    records: list[tuple[int, int, str, bool, int, int]] = []
    current_length = SPACER_LENGTH

    for strand in ("+", "-"):
        for index in range(len(PROFILE_LENGTHS) + 1):
            sequence = coding_orf if strand == "+" else _reverse_complement(coding_orf)
            start = current_length + 1
            stop = current_length + len(sequence)
            profile_index = min(index, len(PROFILE_LENGTHS) - 1)
            records.append(
                (
                    start,
                    stop,
                    strand,
                    index < len(PROFILE_LENGTHS),
                    PROFILE_LENGTHS[profile_index],
                    PROFILE_COUNTS[profile_index],
                )
            )
            pieces.extend((sequence, "N" * SPACER_LENGTH))
            current_length += len(sequence) + SPACER_LENGTH

    genome = "".join(pieces)
    return genome, tuple(FixtureOrf(*record) for record in records)


def _write_bam(path: Path, genome_length: int, records: tuple[FixtureOrf, ...]) -> None:
    """Materialize reads whose A-sites encode the synthetic S-curve."""

    path.parent.mkdir(parents=True)
    alignments = []
    for record in records:
        for position in range(record.signal_length):
            if record.strand == "+":
                asite = record.start - 1 + position
                reference_start = asite - 17
            else:
                asite = record.stop - 1 - position
                reference_start = asite - 12
            alignments.extend(
                (reference_start, record.strand == "-")
                for _ in range(record.signal_count)
            )

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": genome_length}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for index, (reference_start, reverse) in enumerate(sorted(alignments)):
            read = pysam.AlignedSegment()
            read.query_name = f"synthetic_{index:06d}"
            read.query_sequence = "A" * 30
            read.query_qualities = pysam.qualitystring_to_array("I" * 30)
            read.flag = 16 if reverse else 0
            read.reference_id = 0
            read.reference_start = reference_start
            read.mapping_quality = 255
            read.cigartuples = [(0, 30)]
            read.set_tag("NH", 1)
            output.write(read)
    pysam.index(str(path))


def _write_inputs(workdir: Path) -> tuple[Path, tuple[FixtureOrf, ...]]:
    """Create a complete predictor fixture without mapping real reads."""

    inputs = workdir / "input data"
    inputs.mkdir(parents=True)

    genome_sequence, records = _build_genome()
    genome = inputs / "tiny synthetic genome.fa"
    genome.write_text(f">{CONTIG} DeepRibo container smoke\n{genome_sequence}\n")

    annotation = inputs / "tiny synthetic annotation.gff"
    gff_rows = ["##gff-version 3\n"]
    for index, record in enumerate(records):
        if not record.annotated:
            continue
        gff_rows.append(
            f"{CONTIG}\tfixture\tCDS\t{record.start}\t{record.stop}\t.\t"
            f"{record.strand}\t0\tID=fixture_cds_{index:03d};\n"
        )
    annotation.write_text("".join(gff_rows))

    fastq = inputs / "synthetic reads.fastq.gz"
    fastq.write_bytes(
        gzip.compress(
            b"@synthetic\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
            mtime=0,
        )
    )
    samples = inputs / "samples.tsv"
    samples.write_text(
        "method\tcondition\treplicate\tfastqFile\tfastqFile2\n"
        f"RIBO\t{CONDITION}\t{REPLICATE}\t{fastq}\t\n"
    )

    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["biologySettings"].update(
        {"genome": str(genome), "annotation": str(annotation), "samples": str(samples)}
    )
    config["predictionSettings"]["deepribo"] = "on"
    config["workflowSettings"]["stages"] = ["predictions"]
    config_path = inputs / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    _write_bam(
        workdir / "maplink" / f"RIBO-{CONDITION}-{REPLICATE}.bam",
        len(genome_sequence),
        records,
    )
    return config_path, records


def _stream(command: list[str], environment: dict[str, str], log: Path) -> str:
    """Run a command while retaining its combined output for assertions."""

    rendered: list[str] = []
    shown = "+ " + " ".join(command)
    print(shown, flush=True)
    with log.open("a", encoding="utf-8") as log_handle:
        log_handle.write(shown + "\n")
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=environment,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log_handle.write(line)
            log_handle.flush()
            rendered.append(line)
        returncode = process.wait()
    output = "".join(rendered)
    if returncode != 0:
        raise RuntimeError(f"command failed with exit status {returncode}")
    return output


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _bedgraph_rows(path: Path) -> list[tuple[str, int, int, int]]:
    return [
        (fields[0], int(fields[1]), int(fields[2]), int(fields[3]))
        for line in path.read_text().splitlines()
        if line
        for fields in [line.split("\t")]
    ]


def _validate_tracks(workdir: Path, records: tuple[FixtureOrf, ...]) -> None:
    directory = workdir / "coverage_deepribo"
    stem = f"{CONDITION}-{REPLICATE}"
    for suffix, strand in (("fwd", "+"), ("rev", "-")):
        actual = _bedgraph_rows(directory / f"{stem}_asite_{suffix}.bedgraph")
        expected = []
        for record in records:
            if record.strand != strand:
                continue
            positions = (
                range(record.start - 1, record.start - 1 + record.signal_length)
                if strand == "+"
                else range(record.stop - record.signal_length, record.stop)
            )
            expected.extend(
                (CONTIG, position, position + 1, record.signal_count)
                for position in positions
            )
        if actual != sorted(expected):
            raise AssertionError(f"production A-site geometry differs on {strand}")

        coverage_rows = _bedgraph_rows(directory / f"{stem}_cov_{suffix}.bedgraph")
        if not coverage_rows:
            raise AssertionError(f"production read coverage is empty on {strand}")
        previous_end = -1
        for contig, start, stop, count in coverage_rows:
            if contig != CONTIG or start < previous_end or stop <= start or count <= 0:
                raise AssertionError(f"malformed production coverage row on {strand}")
            previous_end = stop


def _validate_png(path: Path) -> None:
    data = path.read_bytes()
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise AssertionError(f"{path} lacks the PNG signature")
    offset = 8
    chunk_types = []
    while offset < len(data):
        if offset + 12 > len(data):
            raise AssertionError(f"{path} has a truncated PNG chunk")
        length = int.from_bytes(data[offset : offset + 4], "big")
        chunk_type = data[offset + 4 : offset + 8]
        end = offset + 12 + length
        if end > len(data):
            raise AssertionError(f"{path} has a truncated PNG payload")
        payload = data[offset + 8 : offset + 8 + length]
        checksum = int.from_bytes(data[offset + 8 + length : end], "big")
        if zlib.crc32(chunk_type + payload) & 0xFFFFFFFF != checksum:
            raise AssertionError(f"{path} has a corrupt PNG chunk")
        chunk_types.append(chunk_type)
        offset = end
    if offset != len(data) or not chunk_types or chunk_types[0] != b"IHDR":
        raise AssertionError(f"{path} is not a structurally complete PNG")
    if chunk_types[-1] != b"IEND" or chunk_types.count(b"IEND") != 1:
        raise AssertionError(f"{path} lacks one final PNG IEND chunk")


def _validate_parser(
    workdir: Path, records: tuple[FixtureOrf, ...]
) -> list[dict[str, str]]:
    parsed = workdir / "deepribo" / "parsed" / f"{CONDITION}-{REPLICATE}"
    rows = _read_csv(parsed / "data_list.csv")
    if not rows or tuple(rows[0]) != PARSED_COLUMNS or len(rows) != len(records):
        raise AssertionError("the pinned parser produced an unexpected table shape")

    expected = {record.parser_locus: record for record in records}
    if {row["locus"] for row in rows} != set(expected):
        raise AssertionError("the parser did not recover the exact synthetic ORFs")
    for row in rows:
        record = expected[row["locus"]]
        if row["strand"] != record.strand or int(row["label"]) != record.annotated:
            raise AssertionError(
                f"incorrect strand/annotation label for {row['locus']}"
            )
        for column, wanted in (
            ("coverage", record.read_coverage),
            ("coverage_elo", record.asite_coverage),
            ("rpk", record.read_occupancy),
            ("rpk_elo", record.asite_occupancy),
        ):
            if not math.isclose(float(row[column]), wanted, rel_tol=0, abs_tol=1e-12):
                raise AssertionError(f"incorrect {column} for {row['locus']}")
        for column in ("filename", "filename_counts"):
            artifact = workdir / "deepribo" / "parsed" / row[column]
            if not artifact.is_file() or artifact.stat().st_size == 0:
                raise AssertionError(f"missing parser tensor: {artifact}")
    return rows


def _validate_cutoffs(workdir: Path) -> tuple[float, float]:
    directory = workdir / "deepribo" / "cutoffs" / f"{CONDITION}-{REPLICATE}"
    fields = (directory / "parameters.txt").read_text().strip().split(",")
    if len(fields) != 2:
        raise AssertionError(f"unexpected DeepRibo cutoff fields: {fields!r}")
    rpkm, coverage = map(float, fields)
    if not math.isfinite(rpkm) or rpkm <= 0:
        raise AssertionError(f"invalid occupancy cutoff: {rpkm}")
    if not math.isfinite(coverage) or not 0 <= coverage <= 0.60:
        raise AssertionError(f"invalid coverage cutoff: {coverage}")
    _validate_png(directory / "s_curve.png")

    receipt = (directory / ".complete").read_text().splitlines()
    expected_receipt = [
        "HRIBO DeepRibo cutoff publication v1",
        str((directory / "parameters.txt").resolve()),
        str((directory / "s_curve.png").resolve()),
    ]
    if receipt != expected_receipt:
        raise AssertionError("cutoff receipt does not identify the published pair")
    return rpkm, coverage


def _validate_predictions(
    workdir: Path,
    parsed_rows: list[dict[str, str]],
    cutoffs: tuple[float, float],
    records: tuple[FixtureOrf, ...],
) -> None:
    path = workdir / "deepribo" / f"{CONDITION}-{REPLICATE}" / "predictions.csv"
    rows = _read_csv(path)
    if not rows or not set(PREDICTION_COLUMNS).issubset(rows[0]):
        raise AssertionError("DeepRibo predictions have an unexpected schema")

    rpkm, coverage = cutoffs
    expected = {
        row["locus"]
        for row in parsed_rows
        if float(row["rpk_elo"]) >= rpkm
        and float(row["coverage_elo"]) >= coverage
        and abs(int(row["start_site"]) - int(row["stop_site"])) > 30
    }
    if len({row["locus"] for row in rows}) != len(rows):
        raise AssertionError("DeepRibo predictions contain duplicate loci")
    if {row["locus"] for row in rows} != expected:
        raise AssertionError("model inference did not apply the published cutoffs")
    if {row["strand"] for row in rows} != {"+", "-"}:
        raise AssertionError("the model did not retain candidates on both strands")
    if {row["label"] for row in rows} != {"True", "False"}:
        raise AssertionError("the fixture did not retain annotated and novel ORFs")

    logits = [float(row["pred"]) for row in rows]
    if not all(math.isfinite(value) for value in logits):
        raise AssertionError("DeepRibo emitted a non-finite raw model score")
    if max(logits) - min(logits) <= 1e-8:
        raise AssertionError("DeepRibo emitted a constant raw model score")
    expected_ranks = set(range(len(rows)))
    for column in ("pred_rank", "SS_pred_rank"):
        ranks = {int(float(row[column])) for row in rows}
        if ranks != expected_ranks:
            raise AssertionError(f"DeepRibo {column} is not a complete rank set")
    if any(row["SS"] != "True" for row in rows):
        raise AssertionError("unique-stop fixture candidates were not marked SS")

    full_coordinates = {record.parser_locus: record.identifier for record in records}
    prediction_by_identifier = {full_coordinates[row["locus"]]: row for row in rows}
    gff = workdir / "deepribo" / f"{CONDITION}-{REPLICATE}.deepribo.gff"
    gff_rows = [
        line.split("\t")
        for line in gff.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    expected_identifiers = {full_coordinates[row["locus"]] for row in rows}
    actual_identifiers = set()
    if len(gff_rows) != len(expected_identifiers):
        raise AssertionError("DeepRibo GFF has missing or duplicate rows")
    for row in gff_rows:
        if len(row) != 9 or row[1:3] != ["deepribo", "CDS"] or row[7] != "0":
            raise AssertionError(f"malformed DeepRibo GFF row: {row!r}")
        attribute_fields = row[8].rstrip(";").split(";")
        if any(field.count("=") != 1 for field in attribute_fields):
            raise AssertionError(f"malformed DeepRibo GFF attributes: {row[8]!r}")
        attribute_pairs = [field.split("=", 1) for field in attribute_fields]
        if len({key for key, _ in attribute_pairs}) != len(attribute_pairs):
            raise AssertionError(f"duplicate DeepRibo GFF attribute: {row[8]!r}")
        attributes = dict(attribute_pairs)
        identifier = attributes["ID"]
        if identifier in actual_identifiers:
            raise AssertionError(f"duplicate DeepRibo GFF ID: {identifier}")
        actual_identifiers.add(identifier)
        prediction = prediction_by_identifier.get(identifier)
        if prediction is None:
            raise AssertionError(f"unexpected DeepRibo GFF ID: {identifier}")
        if (
            attributes.get("condition") != CONDITION
            or attributes.get("replicate") != REPLICATE
            or attributes.get("method") != "deepribo"
            or attributes.get("deepribo_distance")
            != str(int(float(prediction["dist"])))
        ):
            raise AssertionError(f"incorrect DeepRibo GFF provenance: {attributes!r}")
        if identifier != f"{row[0]}:{row[3]}-{row[4]}:{row[6]}":
            raise AssertionError(
                f"DeepRibo GFF ID disagrees with its coordinates: {row!r}"
            )
        score = float(row[5])
        if not math.isfinite(score):
            raise AssertionError("DeepRibo GFF contains a non-finite score")
        if not math.isclose(score, float(prediction["pred"]), rel_tol=0, abs_tol=1e-12):
            raise AssertionError("DeepRibo GFF score disagrees with its prediction")
    if actual_identifiers != expected_identifiers:
        raise AssertionError("DeepRibo GFF did not preserve inferred ORF coordinates")
    validator = shutil.which("gt")
    if validator is None:
        raise RuntimeError("GenomeTools is required to validate the DeepRibo GFF")
    subprocess.run(
        [validator, "gff3validator", str(gff)],
        check=True,
        capture_output=True,
        text=True,
    )


def _validate_outputs(workdir: Path, records: tuple[FixtureOrf, ...]) -> None:
    model = workdir / "deepribo" / "DeepRibo_model_v1.pt"
    if model.stat().st_size != 8_292_882:
        raise AssertionError("the fetched DeepRibo model has the wrong size")
    if hashlib.sha256(model.read_bytes()).hexdigest() != (
        "e3742ed7666f07eb72a28c35c70f2c62ac38f1d201e082b69b4f8cf35fbe4aa3"
    ):
        raise AssertionError("the fetched DeepRibo model has the wrong SHA-256")

    engine = workdir / "deepribo" / "s_curve_cutoff_estimation.R"
    engine_text = engine.read_text()
    for marker in (
        "if (fit_idx < 5L)",
        "DeepRibo S-curve produced a non-finite bend",
        "Log mean A-site occupancy per nucleotide",
    ):
        if marker not in engine_text:
            raise AssertionError(f"patched S-curve engine lacks marker: {marker}")

    _validate_tracks(workdir, records)
    parsed_rows = _validate_parser(workdir, records)
    cutoffs = _validate_cutoffs(workdir)
    _validate_predictions(workdir, parsed_rows, cutoffs, records)


def _snapshot(workdir: Path) -> dict[str, tuple[str, int]]:
    roots = (
        workdir / "genomes" / "genome.fa",
        workdir / "annotation" / "annotation.gff",
        workdir / "annotation" / "annotation_processed.gff",
        workdir / "maplink" / f"RIBO-{CONDITION}-{REPLICATE}.bam",
        workdir / "maplink" / f"RIBO-{CONDITION}-{REPLICATE}.bam.bai",
        workdir / "coverage_deepribo",
        workdir / "deepribo" / "DeepRibo_model_v1.pt",
        workdir / "deepribo" / "s_curve_cutoff_estimation.R",
        workdir / "deepribo" / "parsed" / f"{CONDITION}-{REPLICATE}",
        workdir / "deepribo" / "cutoffs" / f"{CONDITION}-{REPLICATE}",
        workdir / "deepribo" / f"{CONDITION}-{REPLICATE}" / "predictions.csv",
        workdir / "deepribo" / f"{CONDITION}-{REPLICATE}.deepribo.gff",
    )
    paths = []
    for root in roots:
        paths.extend(sorted(root.rglob("*")) if root.is_dir() else [root])
    return {
        str(path.relative_to(workdir)): (
            hashlib.sha256(path.read_bytes()).hexdigest(),
            path.stat().st_mtime_ns,
        )
        for path in paths
        if path.is_file()
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
        "--conda-prefix", type=Path, default=REPO / ".snakemake" / "ci-deepribo-conda"
    )
    parser.add_argument(
        "--apptainer-prefix",
        type=Path,
        default=REPO / ".snakemake" / "ci-deepribo-apptainer",
    )
    parser.add_argument("--snakemake", default="snakemake")
    return parser


def main() -> None:
    args = _parser().parse_args()
    conda_prefix = args.conda_prefix.expanduser().resolve()
    apptainer_prefix = args.apptainer_prefix.expanduser().resolve()
    active_prefix = os.environ.get("CONDA_PREFIX")
    if active_prefix is None:
        raise RuntimeError("activate the HRIBO development environment first")
    interpreter_bin = Path(sys.executable).resolve().parent
    try:
        interpreter_bin.relative_to(Path(active_prefix).resolve())
    except ValueError as error:
        raise RuntimeError(
            "the smoke script's Python must belong to the active conda environment"
        ) from error

    requested_snakemake = Path(args.snakemake).expanduser()
    if requested_snakemake.is_absolute():
        snakemake = shutil.which(str(requested_snakemake))
    else:
        adjacent = interpreter_bin / requested_snakemake
        snakemake = (
            str(adjacent)
            if adjacent.is_file() and os.access(adjacent, os.X_OK)
            else shutil.which(args.snakemake)
        )
        if snakemake is not None:
            try:
                Path(snakemake).resolve().relative_to(Path(active_prefix).resolve())
            except ValueError as error:
                raise RuntimeError(
                    "Snakemake must belong to the active conda environment"
                ) from error
    if snakemake is None:
        raise RuntimeError(f"Snakemake executable is unavailable: {args.snakemake}")
    apptainer = shutil.which("apptainer")
    if apptainer is None:
        raise RuntimeError("Apptainer is unavailable on PATH")

    version = subprocess.run(
        [apptainer, "version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    print(f"Using Apptainer {version}", flush=True)
    workdir = _empty_external_workdir(args.workdir)
    config_path, records = _write_inputs(workdir)
    conda_prefix.mkdir(parents=True, exist_ok=True)
    apptainer_prefix.mkdir(parents=True, exist_ok=True)

    cache_sibling = workdir.parent / f".{workdir.name}-runtime-cache"
    environment = {
        **os.environ,
        "PATH": os.pathsep.join([str(interpreter_bin), os.environ.get("PATH", "")]),
        "XDG_CACHE_HOME": os.environ.get("XDG_CACHE_HOME", str(cache_sibling / "xdg")),
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

    target = f"deepribo/{CONDITION}-{REPLICATE}.deepribo.gff"
    allowed_rules = (
        "retrieveGenome",
        "retrieveAnnotation",
        "checkAnnotation",
        "asiteOccupancy",
        "coverage",
        "deepriboGetModel",
        "prepareDeepRiboSCurveScript",
        "parseDeepRibo",
        "parameterEstimation",
        "predictDeepRibo",
        "deepriboGFF",
    )
    command = [
        snakemake,
        target,
        "--cores",
        "2",
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
        *allowed_rules,
        "--snakefile",
        str(SNAKEFILE),
        "--directory",
        str(workdir),
        "--configfile",
        str(config_path),
    ]

    log = workdir / "smoke-command.log"
    first = _stream(command, environment, log)
    for rule in allowed_rules:
        if f"rule {rule}:" not in first:
            raise AssertionError(f"production rule did not execute: {rule}")
    for helper in (
        "patch_deepribo_scurve.py",
        "deepribo_data_parser.py",
        "run_parameter_estimation.py",
        "parameter_estimation.R",
    ):
        cached_command = any(
            "source-cache" in line and helper in line for line in first.splitlines()
        )
        if f"{helper} (cached)" not in first or not cached_command:
            raise AssertionError(f"container helper was not source-cached: {helper}")
    _validate_outputs(workdir, records)
    before = _snapshot(workdir)

    second = _stream(command, environment, log)
    if "Nothing to be done" not in second:
        raise AssertionError("second Snakemake run was not a no-op")
    if _snapshot(workdir) != before:
        raise AssertionError("a no-op rerun changed a published DeepRibo artifact")

    print("DeepRibo Snakemake + Apptainer smoke passed", flush=True)


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
        print(f"DeepRibo container smoke failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
