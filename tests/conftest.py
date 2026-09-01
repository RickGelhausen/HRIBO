"""Shared fixtures for the HRIBO test suite."""

import random
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parent.parent / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))


CONTIGS = {"NC_000913.3": 6000, "pPlasmid1": 2000}


@pytest.fixture
def genome_file(tmp_path: Path) -> Path:
    """A two-contig nucleotide FASTA whose headers carry a description."""
    rng = random.Random(0)
    path = tmp_path / "genome.fa"
    with path.open("w") as handle:
        for name, length in CONTIGS.items():
            handle.write(f">{name} Escherichia coli complete sequence\n")
            sequence = "".join(rng.choice("ACGT") for _ in range(length))
            for start in range(0, length, 70):
                handle.write(sequence[start : start + 70] + "\n")
    return path


@pytest.fixture
def annotation_file(tmp_path: Path) -> Path:
    """A GFF3 annotation consistent with `genome_file`."""
    lines = ["##gff-version 3"]
    index = 0
    for name, length in CONTIGS.items():
        position = 200
        while position + 400 < length:
            index += 1
            feature = "rRNA" if index % 7 == 0 else ("tRNA" if index % 11 == 0 else "CDS")
            strand = "+" if index % 2 else "-"
            lines.append(
                "\t".join(
                    [
                        name,
                        "RefSeq",
                        feature,
                        str(position),
                        str(position + 299),
                        ".",
                        strand,
                        "0",
                        f"ID=gene{index};locus_tag=b{index:04d}",
                    ]
                )
            )
            position += 450
    path = tmp_path / "annotation.gff"
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def samples(tmp_path: Path):
    """A sample sheet with matched RIBO/RNA libraries and existing fastq files."""
    import pandas as pd

    fastq_dir = tmp_path / "fastq"
    fastq_dir.mkdir()
    rows = []
    for method in ("RIBO", "RNA"):
        for condition in ("A", "B"):
            for replicate in ("1", "2"):
                name = f"{method}-{condition}-{replicate}_R1.fastq.gz"
                import gzip

                with gzip.open(fastq_dir / name, "wt") as handle:
                    handle.write("@r1\nACGTACGTAC\n+\nIIIIIIIIII\n")
                rows.append(
                    {
                        "method": method,
                        "condition": condition,
                        "replicate": replicate,
                        "fastqFile": str(fastq_dir / name),
                        "fastqFile2": None,
                    }
                )
    return pd.DataFrame(rows)


@pytest.fixture
def config(genome_file: Path, annotation_file: Path):
    """A configuration mirroring the shipped template."""
    return {
        "biologySettings": {
            "adapterS3": "",
            "adapterS5": "",
            "adapterP3R1": "",
            "adapterP5R1": "",
            "adapterP3R2": "",
            "adapterP5R2": "",
            "genome": str(genome_file),
            "annotation": str(annotation_file),
            "samples": "config/samples.tsv",
            "alternativeStartCodons": ["GTG", "TTG"],
        },
        "differentialExpressionSettings": {
            "features": ["CDS", "sRNA"],
            "contrasts": [],
            "padjCutoff": 0.05,
            "log2fcCutoff": 1.0,
        },
        "predictionSettings": {"deepribo": "on"},
        "readstatSettings": {"readLengths": "10-80"},
        "metageneSettings": {
            "positionsOutsideORF": 100,
            "positionsInORF": 150,
            "filteringMethods": ["overlap", "length", "rpkm"],
            "neighboringGenesDistance": 50,
            "rpkmThreshold": 10.0,
            "lengthCutoff": 50,
            "mappingMethods": ["fiveprime", "threeprime"],
            "readLengths": "25-34",
            "normalizationMethods": ["raw", "cpm", "window"],
            "outputFormats": ["interactive", "svg"],
            "includePlotlyJS": "integrated",
            "colorList": [],
        },
        "workflowSettings": {"stages": "full"},
    }
