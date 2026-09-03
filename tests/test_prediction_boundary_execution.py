"""Execute the prediction branch around, but not through, its predictor engines.

The raw Reparation and DeepRibo results are synthetic fixtures in the programs'
recorded output formats. Everything downstream of those files, plus the
DeepRibo BAM-to-bedgraph input preparation, runs through the production
Snakemake rules with their real Conda tools. This is deliberately not described
as scientific predictor execution: the Reparation engine, DataParser, S-curve
estimation, and neural-network model remain covered by a separate future
container smoke test.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pysam
import yaml


REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"
CONTIG = "fixture_contig"
CONTIG_LENGTH = 600
GFF3_HEADER = "##gff-version 3\n"

KNOWN = (100, 189, "+")
NOVEL = (250, 339, "+")
NEGATIVE = (400, 489, "-")

ALLOWED_RULES = (
    "all",
    "retrieveGenome",
    "retrieveAnnotation",
    "checkAnnotation",
    "asiteOccupancy",
    "coverage",
    "reparationGFF",
    "concatReparation",
    "mergeConditions",
    "mergeAll",
    "filterAll",
    "reannotatedORFs",
    "deepriboGFF",
    "concatDeepRibo",
    "allDeepRibo",
    "filterDeepRibo",
    "readCounts",
    "mapReadsToAnnotation",
    "mappedReadSummary",
    "createExcelSummary",
    "createExcelSummaryDeepRibo",
    "updatedAnnotation",
)

EXCLUDED_ENGINE_RULES = (
    "trim_single",
    "map",
    "uniprotDBRetrieve",
    "reparation",
    "deepriboGetModel",
    "parseDeepRibo",
    "parameterEstimation",
    "predictDeepRibo",
)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _genome_sequence() -> str:
    sequence = list("A" * CONTIG_LENGTH)
    coding_sequence = "ATG" + "GCT" * 28 + "TAA"
    assert len(coding_sequence) == 90
    for start, end, strand in (KNOWN, NOVEL, NEGATIVE):
        planted = coding_sequence if strand == "+" else _reverse_complement(
            coding_sequence
        )
        assert end - start + 1 == len(planted)
        sequence[start - 1 : end] = planted
    return "".join(sequence)


def _write_fastq(path: Path, name: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(f"@{name}\n{'A' * 30}\n+\n{'I' * 30}\n")


def _write_bam(path: Path, reads: list[tuple[int, bool]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for index, (start, reverse) in enumerate(sorted(reads)):
            record = pysam.AlignedSegment()
            record.query_name = f"read-{index}"
            record.query_sequence = "A" * 30
            record.query_qualities = pysam.qualitystring_to_array("I" * 30)
            record.flag = 16 if reverse else 0
            record.reference_id = 0
            record.reference_start = start
            record.mapping_quality = 60
            record.cigartuples = [(0, 30)]
            record.set_tag("NH", 1)
            output.write(record)
    pysam.index(str(path))


REPARATION_COLUMNS = (
    "ORF_locus",
    "strand",
    "length",
    "start_codon",
    "ribo_count",
    "ribo_rpkm",
    "ribo_coverage",
    "SD_score",
    "SD_pos",
    "prob",
    "ORF_type",
    "Reference",
    "Distance_from_aTIS",
)


def _write_reparation_output(path: Path, replicate: str) -> None:
    probabilities = {
        "1": ("0.70", "0.60"),
        "2": ("0.90", "0.80"),
    }
    known_probability, novel_probability = probabilities[replicate]
    rows = [
        (
            f"{CONTIG}:100-186",
            "+",
            "87",
            "ATG",
            "12",
            "1.5",
            "0.9",
            "5.0",
            "-8",
            known_probability,
            "annotated",
            "aTIS",
            "0",
        ),
        (
            f"{CONTIG}:250-336",
            "+",
            "87",
            "ATG",
            "8",
            "1.0",
            "0.8",
            "4.0",
            "-7",
            novel_probability,
            "sORF",
            "novel",
            "-1",
        ),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\t".join(REPARATION_COLUMNS)
        + "\n"
        + "".join("\t".join(row) + "\n" for row in rows)
    )


DEEPRIBO_COLUMNS = (
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
    "pred",
    "pred_rank",
    "SS",
    "dist",
    "SS_pred_rank",
)


def _deepribo_row(
    replicate: str,
    name: str,
    locus: str,
    strand: str,
    prediction: str,
    rank: str,
    distance: str,
    label: str,
) -> tuple[str, ...]:
    return (
        f"A-{replicate}/0/{name}_seq.pt",
        f"A-{replicate}/0/{name}_reads.pt",
        label,
        "False",
        strand,
        "0.5",
        "0.5",
        "1.0",
        "1.0",
        locus.split(":", 1)[1].split("-", 1)[0],
        "ATG",
        locus.rsplit("-", 1)[1],
        "TAA",
        locus,
        "MA",
        "ATGGCTTAA",
        prediction,
        rank,
        "True",
        distance,
        rank,
    )


def _write_deepribo_output(path: Path, replicate: str) -> None:
    scores = {
        "1": ("0.70", "0.60", "-0.20"),
        "2": ("0.90", "0.80", "-0.10"),
    }
    known_score, novel_score, negative_score = scores[replicate]
    rows = [
        _deepribo_row(
            replicate,
            "known",
            f"{CONTIG}:100-187",
            "+",
            known_score,
            "1",
            "0",
            "True",
        ),
        _deepribo_row(
            replicate,
            "novel",
            f"{CONTIG}:250-337",
            "+",
            novel_score,
            "2",
            "-1",
            "False",
        ),
        _deepribo_row(
            replicate,
            "negative",
            f"{CONTIG}:402-489",
            "-",
            negative_score,
            "3",
            "-1",
            "False",
        ),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        ",".join(DEEPRIBO_COLUMNS)
        + "\n"
        + "".join(",".join(row) + "\n" for row in rows)
    )


def _build_project(tmp_path: Path) -> tuple[Path, Path]:
    workdir = tmp_path / "prediction boundary with spaces"
    inputs = workdir / "input data"
    inputs.mkdir(parents=True)

    genome = inputs / "tiny genome.fa"
    genome.write_text(f">{CONTIG} predictor-boundary fixture\n{_genome_sequence()}\n")
    annotation = inputs / "tiny annotation.gff"
    annotation.write_text(
        GFF3_HEADER
        + f"{CONTIG}\tfixture\tgene\t100\t189\t.\t+\t.\t"
        "ID=reference_gene;locus_tag=known;Name=known_gene;\n"
        + f"{CONTIG}\tfixture\tCDS\t100\t189\t.\t+\t0\t"
        "ID=reference_cds;Parent=reference_gene;locus_tag=known;"
        "Name=known_gene;\n"
    )

    sample_rows = []
    bam_reads = {
        "1": [(105, False), (105, False), (260, False), (430, True)],
        "2": [(105, False), (260, False), (260, False), (430, True)],
    }
    for replicate in ("1", "2"):
        fastq = inputs / f"replicate {replicate}.fastq.gz"
        _write_fastq(fastq, f"replicate-{replicate}")
        sample_rows.append(
            f"RIBO\tA\t{replicate}\t{fastq}\t\n"
        )

        bam = workdir / f"maplink/RIBO-A-{replicate}.bam"
        _write_bam(bam, bam_reads[replicate])
        _write_reparation_output(
            workdir / f"reparation/A-{replicate}/Predicted_ORFs.txt",
            replicate,
        )
        _write_deepribo_output(
            workdir / f"deepribo/A-{replicate}/predictions.csv",
            replicate,
        )

    samples = inputs / "samples.tsv"
    samples.write_text(
        "method\tcondition\treplicate\tfastqFile\tfastqFile2\n"
        + "".join(sample_rows)
    )

    config = yaml.safe_load((REPO / "config/config.yaml").read_text())
    config["biologySettings"].update(
        {
            "genome": str(genome),
            "annotation": str(annotation),
            "samples": str(samples),
        }
    )
    config["predictionSettings"]["deepribo"] = "on"
    config["workflowSettings"]["stages"] = ["predictions"]
    config_path = inputs / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return workdir, config_path


def _environment(workdir: Path) -> tuple[dict[str, str], Path]:
    configured_prefix = os.environ.get("HRIBO_TEST_CONDA_PREFIX")
    conda_prefix = Path(configured_prefix) if configured_prefix else REPO / ".snakemake/conda"
    conda_prefix = conda_prefix.resolve()
    environment = {
        **os.environ,
        "PATH": os.pathsep.join(
            [str(Path(sys.executable).parent), os.environ.get("PATH", "")]
        ),
        "XDG_CACHE_HOME": str(workdir / ".cache"),
    }
    return environment, conda_prefix


def _command(
    snakemake_command: list[str], workdir: Path, config_path: Path, conda_prefix: Path
) -> list[str]:
    explicit_deepribo_inputs = [
        "coverage_deepribo/A-1_asite_fwd.bedgraph",
        "coverage_deepribo/A-1_asite_rev.bedgraph",
        "coverage_deepribo/A-1_cov_fwd.bedgraph",
        "coverage_deepribo/A-1_cov_rev.bedgraph",
    ]
    return [
        *snakemake_command,
        "all",
        *explicit_deepribo_inputs,
        "--cores",
        "2",
        "--printshellcmds",
        "--show-failed-logs",
        "--rerun-incomplete",
        "--notemp",
        "--software-deployment-method",
        "conda",
        "--conda-prefix",
        str(conda_prefix),
        "--allowed-rules",
        *ALLOWED_RULES,
        "--snakefile",
        str(SNAKEFILE),
        "--directory",
        str(workdir),
        "--configfile",
        str(config_path),
    ]


def _run(command: list[str], environment: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        capture_output=True,
        text=True,
        env=environment,
        timeout=1200,
    )


def _gff_rows(path: Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]


def _attributes(row: list[str]) -> dict[str, str]:
    return {
        key: value
        for field in row[8].rstrip(";").split(";")
        if "=" in field
        for key, value in [field.split("=", 1)]
    }


def _coordinate(row: list[str]) -> tuple[int, int, str]:
    return int(row[3]), int(row[4]), row[6]


def _rows_by_coordinate(path: Path) -> dict[tuple[int, int, str], list[str]]:
    return {_coordinate(row): row for row in _gff_rows(path)}


def _bedgraph_rows(path: Path) -> list[tuple[str, int, int, int]]:
    return [
        (fields[0], int(fields[1]), int(fields[2]), int(fields[3]))
        for line in path.read_text().splitlines()
        if line
        for fields in [line.split("\t")]
    ]


def _output_snapshot(workdir: Path) -> dict[str, bytes]:
    relative_paths = (
        "tracks/reparation_annotated.gff",
        "tracks/deepribo_merged.gff",
        "tracks/deepribo_merged_plus.gff",
        "tracks/updated_annotation.gff",
        "readcounts/reparation_annotation.gff",
        "readcounts/deepribo_annotation.gff",
        "auxiliary/predictions_reparation.xlsx",
        "auxiliary/predictions_deepribo.xlsx",
        "coverage_deepribo/A-1_asite_fwd.bedgraph",
        "coverage_deepribo/A-1_asite_rev.bedgraph",
        "coverage_deepribo/A-1_cov_fwd.bedgraph",
        "coverage_deepribo/A-1_cov_rev.bedgraph",
    )
    return {path: (workdir / path).read_bytes() for path in relative_paths}


def test_prediction_boundary_executes_real_postprocessing_without_predictors(
    snakemake_command,
    tmp_path,
):
    workdir, config_path = _build_project(tmp_path)
    environment, conda_prefix = _environment(workdir)
    command = _command(snakemake_command, workdir, config_path, conda_prefix)

    result = _run(command, environment)
    rendered = result.stdout + result.stderr
    assert result.returncode == 0, rendered

    for rule in (
        "asiteOccupancy",
        "coverage",
        "reparationGFF",
        "filterAll",
        "reannotatedORFs",
        "deepriboGFF",
        "filterDeepRibo",
        "readCounts",
        "mapReadsToAnnotation",
        "mappedReadSummary",
        "createExcelSummary",
        "createExcelSummaryDeepRibo",
        "updatedAnnotation",
    ):
        assert f"rule {rule}:" in rendered
    for rule in EXCLUDED_ENGINE_RULES:
        assert f"rule {rule}:" not in rendered

    assert not (workdir / "uniprotDB").exists()
    assert not (workdir / "deepribo/DeepRibo_model_v1.pt").exists()
    assert not (workdir / "deepribo/parsed").exists()
    assert not (workdir / "logs/A-1_reparation.log").exists()
    assert not (workdir / "logs/A-1_predict_deepribo.log").exists()

    assert _bedgraph_rows(
        workdir / "coverage_deepribo/A-1_asite_fwd.bedgraph"
    ) == [
        (CONTIG, 122, 123, 2),
        (CONTIG, 277, 278, 1),
    ]
    assert _bedgraph_rows(
        workdir / "coverage_deepribo/A-1_asite_rev.bedgraph"
    ) == [(CONTIG, 442, 443, 1)]
    assert _bedgraph_rows(workdir / "coverage_deepribo/A-1_cov_fwd.bedgraph") == [
        (CONTIG, 105, 135, 2),
        (CONTIG, 260, 290, 1),
    ]
    assert _bedgraph_rows(workdir / "coverage_deepribo/A-1_cov_rev.bedgraph") == [
        (CONTIG, 430, 460, 1)
    ]

    reparation = _rows_by_coordinate(workdir / "tracks/reparation.gff")
    assert set(reparation) == {KNOWN, NOVEL}
    assert _attributes(reparation[KNOWN])["prob"] == "0.9"
    assert _attributes(reparation[KNOWN])["evidence"] == (
        "reparation-A-1 reparation-A-2"
    )
    assert _attributes(reparation[NOVEL])["prob"] == "0.8"

    deepribo_main = _rows_by_coordinate(workdir / "tracks/deepribo_merged.gff")
    deepribo_plus = _rows_by_coordinate(
        workdir / "tracks/deepribo_merged_plus.gff"
    )
    assert set(deepribo_main) == {KNOWN, NOVEL, NEGATIVE}
    assert set(deepribo_plus) == {KNOWN, NOVEL}
    assert _attributes(deepribo_main[KNOWN])["pred_value"] == "0.9"
    assert _attributes(deepribo_main[KNOWN])["evidence"] == "A-1 A-2"
    assert _attributes(deepribo_main[NEGATIVE])["pred_value"] == "-0.1"

    updated_rows = _gff_rows(workdir / "tracks/updated_annotation.gff")
    prediction_coordinates = {
        _coordinate(row) for row in updated_rows if row[1] in {"reparation", "deepribo"}
    }
    assert prediction_coordinates == {KNOWN, NOVEL}
    assert NEGATIVE not in {_coordinate(row) for row in updated_rows}
    assert sum(row[1] == "reparation" for row in updated_rows) == 2
    assert sum(row[1] == "deepribo" for row in updated_rows) == 2
    assert all(
        _attributes(row)["ID"].startswith(f"{row[1]}:")
        for row in updated_rows
        if row[1] in {"reparation", "deepribo"}
    )

    validator = shutil.which("gt")
    assert validator is not None, "GenomeTools is required by the development environment"
    validation = subprocess.run(
        [validator, "gff3validator", str(workdir / "tracks/updated_annotation.gff")],
        capture_output=True,
        text=True,
    )
    assert validation.returncode == 0, validation.stderr

    mapped_reparation = _rows_by_coordinate(
        workdir / "readcounts/reparation_annotation.gff"
    )
    mapped_deepribo = _rows_by_coordinate(
        workdir / "readcounts/deepribo_annotation.gff"
    )
    assert mapped_reparation[KNOWN][9:] == ["2", "1"]
    assert mapped_reparation[NOVEL][9:] == ["1", "2"]
    assert mapped_deepribo[KNOWN][9:] == ["2", "1"]
    assert mapped_deepribo[NOVEL][9:] == ["1", "2"]
    assert mapped_deepribo[NEGATIVE][9:] == ["1", "1"]

    reparation_workbook = pd.read_excel(
        workdir / "auxiliary/predictions_reparation.xlsx",
        sheet_name="CDS",
        engine="openpyxl",
    )
    deepribo_workbook = pd.read_excel(
        workdir / "auxiliary/predictions_deepribo.xlsx",
        sheet_name="CDS",
        engine="openpyxl",
    )
    assert len(reparation_workbook) == 2
    assert len(deepribo_workbook) == 3
    for workbook in (reparation_workbook, deepribo_workbook):
        assert {"RIBO-A-1_rpkm", "RIBO-A-2_rpkm", "Evidence"} <= set(
            workbook.columns
        )

    known_reparation = reparation_workbook.loc[
        reparation_workbook["Start"] == KNOWN[0]
    ].iloc[0]
    assert known_reparation["Reparation_probability"] == 0.9
    assert known_reparation["RIBO-A-1_rpkm"] == 5_555_555.56
    assert known_reparation["RIBO-A-2_rpkm"] == 2_777_777.78
    negative_deepribo = deepribo_workbook.loc[
        deepribo_workbook["Start"] == NEGATIVE[0]
    ].iloc[0]
    assert negative_deepribo["Deepribo_score"] == -0.1
    assert negative_deepribo["RIBO-A-1_rpkm"] == 2_777_777.78
    assert negative_deepribo["RIBO-A-2_rpkm"] == 2_777_777.78

    before = _output_snapshot(workdir)
    rerun = _run(command, environment)
    rerun_rendered = rerun.stdout + rerun.stderr
    assert rerun.returncode == 0, rerun_rendered
    assert "Nothing to be done" in rerun_rendered
    assert _output_snapshot(workdir) == before
