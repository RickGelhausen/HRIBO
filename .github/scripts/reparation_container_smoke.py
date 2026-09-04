#!/usr/bin/env python3
"""Exercise HRIBO's production Reparation core from preseeded synthetic inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import zlib
from collections import Counter
from pathlib import Path

import pysam


REPO = Path(__file__).resolve().parents[2]
REPARATION_RULE = REPO / "workflow" / "rules" / "reparation.smk"
CONTAINER_REGISTRY = "quay.io/biocontainers/reparation_blast"
CONTAINER_DIGEST = (
    "6852b3b69b532039a5d674115b9cbfd2953f6479cf187bcc769e2d899fcbc288"
)
CONTAINER = f"docker://{CONTAINER_REGISTRY}@sha256:{CONTAINER_DIGEST}"
GENOME_SHA256 = "7c0ea8122a0f5a05cd8357a8693aaba6cd54375362078cdf856535ad01af5993"
PROTEIN_DB_SHA256 = (
    "c06d9388147d5da67a367d787afeb248bd9fc97ee487433999ab1468669aa3b9"
)
READ_LENGTH = 30
FALLBACK_OFFSET = 13
EXPECTED_PRODIGAL_CALLS = 288
EXPECTED_POSITIVES = 286
EXPECTED_NEGATIVE_CANDIDATES = 48
EXPECTED_READS = 277_764
EXPECTED_PREDICTIONS = 321
ORF_HEADER = [
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
]
CODON_TABLE = {
    "GCC": "A",
    "GAA": "E",
    "GGC": "G",
    "AAG": "K",
    "CCG": "P",
    "TCC": "S",
    "CGC": "R",
    "AAC": "N",
    "GAC": "D",
}
BODY_CODONS = tuple(CODON_TABLE)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _assert_digest(path: Path, expected: str) -> None:
    observed = _sha256(path)
    if observed != expected:
        raise AssertionError(
            f"fixture drift in {path}: expected SHA-256 {expected}, got {observed}"
        )


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _coding_body(index: int, codons: int) -> str:
    sequence = "".join(
        BODY_CODONS[(position + index) % len(BODY_CODONS)]
        for position in range(codons)
    )
    forbidden_codons = ("ATG", "GTG", "TTG", "CTG", "TAA", "TAG", "TGA")
    for forbidden in forbidden_codons:
        if forbidden in sequence:
            raise AssertionError(f"coding body unexpectedly contains {forbidden}")
    return sequence


def _write_genome(workdir: Path) -> Path:
    """Create the stop-bounded two-strand genome used by the exact-image probe."""

    padding = "TAATAATAATAATAATAATAATAATAATAA" + "ACACACAC" + "GGAGGAAAAAA"
    blocks = ["C" * 120]
    for index in range(240):
        gene = "ATG" + _coding_body(index, 145 + index % 9) + "TAA"
        block = padding + gene + "TAATAATAATAATAA"
        blocks.append(block if index % 2 == 0 else _reverse_complement(block))
    for index in range(48):
        gene = "CTG" + _coding_body(500 + index, 155 + index % 7) + "TAA"
        block = padding + gene + "TAATAATAATAATAA"
        blocks.append(block if index % 2 == 0 else _reverse_complement(block))
    blocks.append("G" * 120)

    sequence = "".join(blocks)
    genome = workdir / "genomes" / "genome.fa"
    genome.parent.mkdir(parents=True)
    wrapped = "\n".join(
        sequence[index : index + 80] for index in range(0, len(sequence), 80)
    )
    genome.write_text(f">fixture\n{wrapped}\n")
    _assert_digest(genome, GENOME_SHA256)
    return genome


def _read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text().splitlines()
        if not line.startswith(">")
    )


def _translate(sequence: str) -> str:
    amino_acids = (
        "FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG"
    )
    codons = [
        first + second + third
        for first in "TCAG"
        for second in "TCAG"
        for third in "TCAG"
    ]
    table = dict(zip(codons, amino_acids, strict=True))
    return "".join(
        table[sequence[index : index + 3]] for index in range(0, len(sequence), 3)
    )


def _parse_prodigal(path: Path) -> list[tuple[int, int, str]]:
    calls = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise AssertionError(f"malformed Prodigal GFF row: {line!r}")
        calls.append((int(fields[3]), int(fields[4]), fields[6]))
    return calls


def _orf_coordinates(call: tuple[int, int, str]) -> tuple[int, int, str]:
    start, end, strand = call
    return (start, end - 3, strand) if strand == "+" else (start + 3, end, strand)


def _find_start_orfs(
    genome: str, start_codon: str
) -> list[tuple[int, int, str, str]]:
    result = []
    for sequence, strand in ((genome, "+"), (_reverse_complement(genome), "-")):
        offset = 0
        while True:
            start = sequence.find(start_codon, offset)
            if start < 0:
                break
            offset = start + 1
            stop = None
            for cursor in range(start, len(sequence) - 2, 3):
                if sequence[cursor : cursor + 3] in {"TAA", "TAG", "TGA"}:
                    stop = cursor
                    break
            if stop is None or stop - start < 30:
                continue
            if strand == "+":
                genomic_start, genomic_end = start + 1, stop
            else:
                genomic_start = len(genome) - stop + 1
                genomic_end = len(genome) - start
            result.append((genomic_start, genomic_end, strand, start_codon))
    return sorted(set(result))


def _sample_positions(start: int, end: int, count: int) -> list[int]:
    length = end - start + 1
    if count >= length:
        return list(range(start, end + 1))
    return sorted(
        {start + (index * (length - 1)) // (count - 1) for index in range(count)}
    )


def _write_bam(
    path: Path,
    genome_length: int,
    positives: list[tuple[int, int, str, str]],
    negatives: list[tuple[int, int, str, str]],
) -> int:
    """Write the deterministic logistic coverage fixture as an indexed BAM."""

    records: list[tuple[int, int, int]] = []
    read_index = 0
    covered = positives + negatives
    for index, (start, end, strand, _) in enumerate(covered):
        if index < len(positives):
            x_value = -1.15 + 2.65 * ((index % 60) / 59)
        else:
            x_value = 0.65 + 0.75 * ((index - len(positives)) % 12) / 11
        coverage = 0.10 + 0.84 / (1.0 + math.exp(-(x_value - 0.05) / 0.42))
        density = math.exp(x_value)
        length = end - start + 1
        occupied = max(4, min(length, round(coverage * length)))
        total_reads = max(occupied, round(density * length))
        positions = _sample_positions(start, end, occupied)
        repeats, extra = divmod(total_reads, len(positions))
        for position_index, occupancy_position in enumerate(positions):
            copies = repeats + (position_index < extra)
            if strand == "+":
                alignment_start = occupancy_position - FALLBACK_OFFSET
                flag = 0
            else:
                alignment_start = occupancy_position - (
                    READ_LENGTH - 1 - FALLBACK_OFFSET
                )
                flag = 16
            if (
                alignment_start < 1
                or alignment_start + READ_LENGTH - 1 > genome_length
            ):
                continue
            for _ in range(copies):
                read_index += 1
                records.append((alignment_start, read_index, flag))

    records.sort()
    path.parent.mkdir(parents=True)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "fixture", "LN": genome_length}],
    }
    qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for alignment_start, identifier, flag in records:
            record = pysam.AlignedSegment()
            record.query_name = f"read_{identifier:09d}"
            record.flag = flag
            record.reference_id = 0
            record.reference_start = alignment_start - 1
            record.mapping_quality = 60
            record.cigar = ((0, READ_LENGTH),)
            record.query_sequence = "A" * READ_LENGTH
            record.query_qualities = qualities
            handle.write(record)
    pysam.index(str(path))
    return len(records)


def _finish_fixture(workdir: Path, prodigal_gff: Path) -> None:
    """Derive the tiny homolog DB, checked annotation, and BAM from Prodigal."""

    genome_path = workdir / "genomes" / "genome.fa"
    genome = _read_fasta(genome_path)
    calls = _parse_prodigal(prodigal_gff)
    if len(calls) != EXPECTED_PRODIGAL_CALLS:
        raise AssertionError(
            f"pinned image produced {len(calls)} Prodigal calls, "
            f"expected {EXPECTED_PRODIGAL_CALLS}"
        )

    positives: list[tuple[int, int, str, str]] = []
    proteins: dict[str, str] = {}
    for start, end, strand in map(_orf_coordinates, calls):
        dna = genome[start - 1 : end]
        if strand == "-":
            dna = _reverse_complement(dna)
        start_codon = dna[:3]
        if start_codon not in {"ATG", "GTG", "TTG"}:
            continue
        positives.append((start, end, strand, start_codon))
        protein = _translate(dna)
        proteins.setdefault(protein, f"protein_{len(proteins) + 1:04d}")

    negatives = _find_start_orfs(genome, "CTG")
    if len(positives) != EXPECTED_POSITIVES:
        raise AssertionError(
            f"fixture has {len(positives)} positive calls, expected {EXPECTED_POSITIVES}"
        )
    if len(negatives) != EXPECTED_NEGATIVE_CANDIDATES:
        raise AssertionError(
            f"fixture has {len(negatives)} CTG candidates, "
            f"expected {EXPECTED_NEGATIVE_CANDIDATES}"
        )

    database = workdir / "uniprotDB" / "uniprot_sprot.fasta"
    database.parent.mkdir(parents=True)
    with database.open("w") as handle:
        for protein, identifier in proteins.items():
            handle.write(f">{identifier}\n{protein}\n")
    _assert_digest(database, PROTEIN_DB_SHA256)

    annotation = workdir / "annotation" / "annotation_processed.gff"
    annotation.parent.mkdir(parents=True)
    with annotation.open("w") as handle:
        handle.write("##gff-version 3\n")
        for index, (start, end, strand, _) in enumerate(positives, 1):
            if strand == "+":
                annotated_start, annotated_end = start, end + 3
            else:
                annotated_start, annotated_end = start - 3, end
            identifier = f"gene_{index:04d}"
            attributes = (
                f"ID=cds_{index:04d};locus_tag={identifier};"
                f"Name={identifier};gene_biotype=protein_coding"
            )
            handle.write(
                "\t".join(
                    (
                        "fixture",
                        "fixture",
                        "CDS",
                        str(annotated_start),
                        str(annotated_end),
                        ".",
                        strand,
                        "0",
                        attributes,
                    )
                )
                + "\n"
            )

    reads = _write_bam(
        workdir / "maplink" / "RIBO-A-1.bam",
        len(genome),
        positives,
        negatives,
    )
    if reads != EXPECTED_READS:
        raise AssertionError(f"fixture has {reads} reads, expected {EXPECTED_READS}")
    print(
        "Fixture inputs ready: "
        f"{len(calls)} Prodigal calls, {len(positives)} positives, "
        f"{len(negatives)} CTG candidates, {len(proteins)} unique homologs, "
        f"{reads} alignments",
        flush=True,
    )


def _stream(
    command: list[str], environment: dict[str, str], log_path: Path
) -> str:
    """Stream and durably retain combined output from one Snakemake phase."""

    print("+ " + " ".join(command), flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_handle:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=environment,
        )
        assert process.stdout is not None
        rendered = []
        for line in process.stdout:
            print(line, end="", flush=True)
            log_handle.write(line)
            log_handle.flush()
            rendered.append(line)
        returncode = process.wait()
    output = "".join(rendered)
    if returncode != 0:
        raise RuntimeError(
            f"command failed with exit status {returncode}; see {log_path}"
        )
    return output


def _write_fixture_snakefile(workdir: Path) -> Path:
    path = workdir / "fixture-build" / "ProdigalSnakefile"
    path.parent.mkdir(parents=True)
    path.write_text(
        f'''rule fixtureProdigal:
    input:
        "genomes/genome.fa"
    output:
        "fixture-build/prodigal.gff"
    container:
        {CONTAINER!r}
    threads: 1
    shell:
        "prodigal -i {{input:q}} -f gff -o {{output:q}} -g 11 -q"
'''
    )
    return path


def _write_production_snakefile(workdir: Path) -> Path:
    path = workdir / "fixture-build" / "ReparationSnakefile"
    path.write_text(
        f'''from pathlib import Path
import pandas as pd

SCRIPTS = Path({str(REPO / "workflow" / "scripts")!r})
samples = pd.DataFrame(
    [{{"method": "RIBO", "condition": "A", "replicate": "1"}}]
)
conditions = ["A"]


rule retrieveGenome:
    output:
        "genomes/genome.fa"


rule checkAnnotation:
    output:
        "annotation/annotation_processed.gff"


include: {str(REPARATION_RULE)!r}
'''
    )
    return path


def _deployment_arguments(apptainer_prefix: Path) -> list[str]:
    return [
        "--cores",
        "1",
        "--printshellcmds",
        "--show-failed-logs",
        "--rerun-incomplete",
        "--notemp",
        "--latency-wait",
        "60",
        "--apptainer-prefix",
        str(apptainer_prefix),
    ]


def _validate_production_container_source() -> None:
    source = REPARATION_RULE.read_text()
    expected_registry = f'"docker://{CONTAINER_REGISTRY}@sha256:"'
    expected_digest = f'"{CONTAINER_DIGEST}"'
    if expected_registry not in source or expected_digest not in source:
        raise AssertionError(
            "production Reparation rule no longer names the expected digest-pinned image"
        )


def _validate_cached_image(apptainer: str, apptainer_prefix: Path) -> None:
    cache_name = hashlib.md5(
        CONTAINER.encode(), usedforsecurity=False
    ).hexdigest()
    image = apptainer_prefix / f"{cache_name}.simg"
    if not image.is_file() or image.stat().st_size == 0:
        raise AssertionError(f"Snakemake did not cache the expected image: {image}")
    metadata = json.loads(
        subprocess.run(
            [apptainer, "inspect", "--json", str(image)],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
    )
    labels = metadata["data"]["attributes"]["labels"]
    if labels.get("org.opencontainers.image.base.digest") != (
        f"sha256:{CONTAINER_DIGEST}"
    ):
        raise AssertionError("cached Reparation SIF has an unexpected source digest")
    if labels.get("org.opencontainers.image.base.name") != CONTAINER.removeprefix(
        "docker://"
    ):
        raise AssertionError("cached Reparation SIF has an unexpected source image")


def _png_dimensions(path: Path) -> tuple[int, int]:
    data = path.read_bytes()
    if not data.startswith(b"\x89PNG\r\n\x1a\n") or len(data) < 33:
        raise AssertionError(f"{path} is not a complete PNG")
    length = struct.unpack(">I", data[8:12])[0]
    if length != 13 or data[12:16] != b"IHDR":
        raise AssertionError(f"{path} has no canonical IHDR chunk")
    if zlib.crc32(data[12:29]) & 0xFFFFFFFF != struct.unpack(">I", data[29:33])[0]:
        raise AssertionError(f"{path} has an invalid IHDR checksum")
    return struct.unpack(">II", data[16:24])


def _pdf_pages(path: Path, pdfinfo: str) -> int:
    data = path.read_bytes()
    if not data.startswith(b"%PDF-") or b"%%EOF" not in data[-2048:]:
        raise AssertionError(f"{path} is not a complete PDF")
    result = subprocess.run(
        [pdfinfo, str(path)], capture_output=True, text=True, check=True
    ).stdout
    for line in result.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1])
    raise AssertionError(f"pdfinfo did not report a page count for {path}")


def _parse_attributes(raw_attributes: str) -> dict[str, str]:
    attributes = {}
    for field in raw_attributes.rstrip(";").split(";"):
        if not field:
            continue
        key, separator, value = field.partition("=")
        if not separator or not key or key in attributes:
            raise AssertionError(f"malformed GFF3 attributes: {raw_attributes!r}")
        attributes[key] = value
    return attributes


def _parse_locus(locus: str) -> tuple[str, int, int]:
    contig, separator, interval = locus.rpartition(":")
    match = re.fullmatch(r"(\d+)-(\d+)", interval)
    if not separator or not contig or match is None:
        raise AssertionError(f"malformed Reparation locus: {locus!r}")
    return contig, int(match.group(1)), int(match.group(2))


def _validate_annotation_adapter(workdir: Path) -> None:
    source = workdir / "annotation" / "annotation_processed.gff"
    source_geometry = set()
    for line in source.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        source_geometry.add((fields[0], int(fields[3]), int(fields[4]), fields[6]))

    adapter = workdir / "reparation" / "annotation.gtf"
    adapter_geometry = set()
    transcript_ids = set()
    for line in adapter.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] != "transcript" or fields[7] != ".":
            raise AssertionError(f"malformed Reparation adapter row: {line!r}")
        geometry = (fields[0], int(fields[3]), int(fields[4]), fields[6])
        adapter_geometry.add(geometry)
        attributes = dict(re.findall(r'(\w+) "([^"]*)";', fields[8]))
        expected_keys = {"gene_id", "transcript_id", "gene_name", "gene_biotype"}
        if set(attributes) != expected_keys or attributes["gene_biotype"] != "protein_coding":
            raise AssertionError(f"unexpected Reparation adapter attributes: {fields[8]!r}")
        transcript_ids.add(attributes["transcript_id"])

    if len(source_geometry) != EXPECTED_POSITIVES or adapter_geometry != source_geometry:
        raise AssertionError("production annotation adapter changed fixture geometries")
    if len(transcript_ids) != EXPECTED_POSITIVES:
        raise AssertionError("production annotation adapter emitted duplicate transcript IDs")


def _table_records(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, object]]:
    records = {}
    for row in rows:
        contig, start, stop = _parse_locus(row["ORF_locus"])
        identity = (row["ORF_locus"], row["strand"])
        if identity in records:
            raise AssertionError(f"duplicate Reparation prediction: {identity!r}")
        records[identity] = {
            "contig": contig,
            "start": start,
            "stop": stop,
            "length": int(row["length"]),
            "codon": row["start_codon"],
            "row": row,
        }
    return records


def _validate_bed(result: Path, expected: dict[tuple[str, str], dict[str, object]]) -> None:
    lines = (result / "Predicted_ORFs.bed").read_text().splitlines()
    if not lines or not lines[0].startswith("track type=bed "):
        raise AssertionError("Predicted_ORFs.bed lacks its track header")
    seen = set()
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != 9:
            raise AssertionError(f"malformed Reparation BED row: {line!r}")
        identity = (fields[3], fields[5])
        if identity not in expected or identity in seen:
            raise AssertionError(f"unexpected or duplicate Reparation BED row: {identity!r}")
        seen.add(identity)
        record = expected[identity]
        start = int(record["start"])
        stop = int(record["stop"])
        if fields[5] == "+":
            coordinates = (start - 1, stop + 3, start - 1, stop)
        else:
            coordinates = (start - 4, stop, start - 1, stop)
        observed = tuple(int(fields[index]) for index in (1, 2, 6, 7))
        if fields[0] != record["contig"] or observed != coordinates:
            raise AssertionError(f"BED coordinates disagree with {identity!r}")
    if seen != set(expected):
        raise AssertionError("Predicted_ORFs.bed does not cover the exact table result set")


def _validate_fasta(
    result: Path, expected: dict[tuple[str, str], dict[str, object]]
) -> None:
    fasta_records = []
    header = None
    sequence = []
    for line in (result / "Predicted_ORFs.fasta").read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                fasta_records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        elif line:
            sequence.append(line)
    if header is not None:
        fasta_records.append((header, "".join(sequence)))

    seen = set()
    pattern = re.compile(
        r"generic\|(.+)\|start codon:([ACGT]{3}) strand:([+-]) length:(\d+)"
    )
    for raw_header, protein in fasta_records:
        match = pattern.fullmatch(raw_header)
        if match is None:
            raise AssertionError(f"malformed Reparation FASTA header: {raw_header!r}")
        identity = (match.group(1), match.group(3))
        if identity not in expected or identity in seen:
            raise AssertionError(
                f"unexpected or duplicate Reparation FASTA record: {identity!r}"
            )
        seen.add(identity)
        record = expected[identity]
        if (
            match.group(2) != record["codon"]
            or int(match.group(4)) != record["length"]
            or len(protein) != int(record["length"]) // 3
            or not re.fullmatch(r"[A-Za-z*.-]+", protein)
        ):
            raise AssertionError(f"FASTA metadata disagrees with {identity!r}")
    if seen != set(expected):
        raise AssertionError("Predicted_ORFs.fasta does not cover the exact table result set")


def _validate_offsets(result: Path) -> None:
    rows = []
    found_header = False
    for raw_line in (result / "p_site_offsets.txt").read_text().splitlines():
        if raw_line == "length\tp_offset":
            found_header = True
            continue
        if found_header and raw_line and not raw_line.startswith("#"):
            length, offset = raw_line.split("\t")
            rows.append((length, int(offset)))
    expected = [(str(length), 13) for length in range(22, 41)] + [("default", 13)]
    if rows != expected:
        raise AssertionError(f"unexpected exact-image P-site offsets: {rows!r}")


def _validate_gff(
    workdir: Path,
    table_rows: list[dict[str, str]],
    strands: Counter[str],
) -> Path:
    gff = workdir / "reparation" / "A-1.reparation.gff"
    gff_lines = [line for line in gff.read_text().splitlines() if line]
    if gff_lines[:1] != ["##gff-version 3"]:
        raise AssertionError("Reparation GFF lacks its version directive")
    records = [line.split("\t") for line in gff_lines[1:] if not line.startswith("#")]
    if len(records) != EXPECTED_PREDICTIONS or any(len(row) != 9 for row in records):
        raise AssertionError("Reparation GFF does not contain 321 nine-column records")
    if Counter(row[6] for row in records) != strands:
        raise AssertionError("Reparation GFF strand counts differ from its result table")

    for source, converted in zip(table_rows, records, strict=True):
        contig, start, stop = _parse_locus(source["ORF_locus"])
        strand = source["strand"]
        converted_start = start if strand == "+" else start - 3
        converted_stop = stop + 3 if strand == "+" else stop
        identifier = f"{contig}:{converted_start}-{converted_stop}:{strand}"
        if converted[:8] != [
            contig,
            "reparation",
            "CDS",
            str(converted_start),
            str(converted_stop),
            ".",
            strand,
            "0",
        ]:
            raise AssertionError(f"GFF geometry disagrees with {source['ORF_locus']!r}")
        attributes = _parse_attributes(converted[8])
        expected_attribute_names = {
            "ID",
            "Name",
            "orf_type",
            "length",
            "ribo_count",
            "ribo_rpkm",
            "ribo_coverage",
            "sd_score",
            "sd_pos",
            "prob",
            "reference",
            "distance_from_atis",
            "condition",
            "replicate",
            "method",
        }
        if set(attributes) != expected_attribute_names:
            raise AssertionError(
                f"GFF attribute names disagree with {source['ORF_locus']!r}"
            )
        exact_attributes = {
            "ID": identifier,
            "Name": identifier,
            "orf_type": source["ORF_type"],
            "condition": "A",
            "replicate": "1",
            "method": "reparation",
        }
        if any(attributes.get(key) != value for key, value in exact_attributes.items()):
            raise AssertionError(f"GFF attributes disagree with {source['ORF_locus']!r}")
        if int(attributes["length"]) != int(source["length"]):
            raise AssertionError(f"GFF length disagrees with {source['ORF_locus']!r}")
        numeric_fields = {
            "ribo_count": "ribo_count",
            "ribo_rpkm": "ribo_rpkm",
            "ribo_coverage": "ribo_coverage",
            "sd_score": "SD_score",
            "prob": "prob",
        }
        for attribute, table_field in numeric_fields.items():
            if not math.isclose(
                float(attributes[attribute]),
                float(source[table_field]),
                rel_tol=1e-12,
                abs_tol=1e-12,
            ):
                raise AssertionError(
                    f"GFF {attribute} disagrees with {source['ORF_locus']!r}"
                )
        for attribute, table_field in (
            ("sd_pos", "SD_pos"),
            ("distance_from_atis", "Distance_from_aTIS"),
        ):
            table_value = source[table_field]
            converted_value = attributes[attribute]
            if table_value == "NA":
                matches = converted_value.lower() == "nan"
            else:
                matches = float(converted_value) == float(table_value)
            if not matches:
                raise AssertionError(
                    f"GFF {attribute} disagrees with {source['ORF_locus']!r}"
                )
        reference = attributes["reference"]
        expected_reference = source["Reference"]
        if expected_reference == "NA":
            reference_matches = reference.lower() == "nan"
        else:
            reference_matches = reference == expected_reference
        if not reference_matches:
            raise AssertionError(
                f"GFF reference disagrees with {source['ORF_locus']!r}"
            )
    return gff


def _validate_outputs(
    workdir: Path, input_snapshot: dict[str, str]
) -> None:
    result = workdir / "reparation" / "A-1"
    receipt = result / ".complete"
    if receipt.read_bytes() != b"HRIBO REPARATION publication v1\n":
        raise AssertionError("production completion receipt has unexpected content")

    with (result / "Predicted_ORFs.txt").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ORF_HEADER:
            raise AssertionError(f"unexpected Reparation header: {reader.fieldnames!r}")
        rows = list(reader)
    if len(rows) != EXPECTED_PREDICTIONS:
        raise AssertionError(
            f"Reparation emitted {len(rows)} predictions, expected {EXPECTED_PREDICTIONS}"
        )
    strands = Counter(row["strand"] for row in rows)
    if strands != Counter({"+": 159, "-": 162}):
        raise AssertionError(f"unexpected prediction strands: {dict(strands)!r}")
    classes = Counter(row["ORF_type"] for row in rows)
    expected_classes = Counter(
        {"Annotated": 226, "5' extension": 57, "Intergenic": 36, "Truncation": 2}
    )
    if classes != expected_classes:
        raise AssertionError(f"unexpected prediction classes: {dict(classes)!r}")

    expected_records = _table_records(rows)
    if len(expected_records) != EXPECTED_PREDICTIONS:
        raise AssertionError("Reparation table does not have unique locus/strand rows")
    _validate_annotation_adapter(workdir)
    _validate_bed(result, expected_records)
    _validate_fasta(result, expected_records)
    _validate_offsets(result)
    gff = _validate_gff(workdir, rows, strands)

    genome_tools = shutil.which("gt")
    if genome_tools is None:
        raise RuntimeError("GenomeTools `gt` is required to validate the produced GFF3")
    subprocess.run([genome_tools, "gff3validator", str(gff)], check=True)

    pdfinfo = shutil.which("pdfinfo")
    if pdfinfo is None:
        raise RuntimeError("pdfinfo is required to validate Reparation reports")
    expected_pages = {
        "metagene_profile.pdf": 2,
        "PR_and_ROC_curve.pdf": 2,
        "variable_importance.pdf": 1,
        "S_Curve.pdf": 1,
    }
    for name, pages in expected_pages.items():
        observed = _pdf_pages(result / name, pdfinfo)
        if observed != pages:
            raise AssertionError(f"{name} has {observed} pages, expected {pages}")
    dimensions = _png_dimensions(result / "p_site_offset.png")
    if dimensions != (930, 2405):
        raise AssertionError(
            f"p_site_offset.png is {dimensions[0]}x{dimensions[1]}, expected 930x2405"
        )

    engine_log = (workdir / "logs" / "A-1_reparation.log").read_text()
    normalized_log = " ".join(engine_log.split())
    for marker in (
        "Total number of ORFs in positive set 286",
        "Number of ORFs in negative Set 27",
        "Total number of ORF families predicted 408",
    ):
        # The first message comes from Perl while the latter messages come from
        # R and differ only in harmless whitespace across builds.
        if marker not in normalized_log:
            raise AssertionError(f"Reparation engine log lacks marker: {marker}")

    database = workdir / "uniprotDB" / "uniprot_sprot.fasta"
    sidecars = [path for path in database.parent.iterdir() if path != database]
    if sidecars:
        raise AssertionError(f"source protein DB was mutated: {sidecars!r}")
    _assert_digest(database, PROTEIN_DB_SHA256)
    _assert_digest(workdir / "genomes" / "genome.fa", GENOME_SHA256)

    after = _input_snapshot(workdir)
    if after != input_snapshot:
        changed = sorted(path for path in after if after[path] != input_snapshot.get(path))
        raise AssertionError(f"production execution mutated fixture inputs: {changed!r}")

    reparation_root = workdir / "reparation"
    debris = [
        reparation_root / name
        for name in (
            ".A-1_candidate",
            ".A-1_backup",
            ".A-1_transaction",
            ".A-1_transaction.tmp",
            ".A-1_staging",
            ".A-1_staging.tmp",
        )
        if os.path.lexists(reparation_root / name)
    ]
    for name in (".complete.tmp",):
        path = result / name
        if os.path.lexists(path):
            debris.append(path)
    if debris:
        raise AssertionError(f"Reparation left publication transaction debris: {debris!r}")


def _input_snapshot(workdir: Path) -> dict[str, str]:
    paths = (
        workdir / "genomes" / "genome.fa",
        workdir / "annotation" / "annotation_processed.gff",
        workdir / "uniprotDB" / "uniprot_sprot.fasta",
        workdir / "maplink" / "RIBO-A-1.bam",
        workdir / "maplink" / "RIBO-A-1.bam.bai",
    )
    return {str(path.relative_to(workdir)): _sha256(path) for path in paths}


def _snapshot(workdir: Path) -> dict[str, tuple[str, int]]:
    result = workdir / "reparation" / "A-1"
    paths = [path for path in result.rglob("*") if path.is_file()]
    paths.extend(
        (
            workdir / "reparation" / "annotation.gtf",
            workdir / "reparation" / "A-1.reparation.gff",
        )
    )
    return {
        str(path.relative_to(workdir)): (_sha256(path), path.stat().st_mtime_ns)
        for path in sorted(paths)
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
        default=REPO / ".snakemake" / "ci-reparation-conda",
    )
    parser.add_argument(
        "--apptainer-prefix",
        type=Path,
        default=REPO / ".snakemake" / "ci-reparation-apptainer",
    )
    parser.add_argument("--snakemake", default="snakemake")
    return parser


def _resolve_executable(command: str) -> str:
    expanded = Path(command).expanduser()
    if expanded.parent != Path("."):
        resolved = expanded.resolve()
        if not resolved.is_file() or not os.access(resolved, os.X_OK):
            raise RuntimeError(f"executable is unavailable: {command}")
        return str(resolved)
    sibling = Path(sys.executable).resolve().parent / command
    if sibling.is_file() and os.access(sibling, os.X_OK):
        return str(sibling)
    resolved = shutil.which(command)
    if resolved is None:
        raise RuntimeError(f"executable is unavailable: {command}")
    return resolved


def main() -> None:
    args = _parser().parse_args()
    workdir = _empty_external_workdir(args.workdir)
    conda_prefix = args.conda_prefix.expanduser().resolve()
    apptainer_prefix = args.apptainer_prefix.expanduser().resolve()
    snakemake = _resolve_executable(args.snakemake)
    apptainer = shutil.which("apptainer")
    if apptainer is None:
        raise RuntimeError("Apptainer is unavailable on PATH")
    version = subprocess.run(
        [apptainer, "version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    print(f"Using Apptainer {version}", flush=True)

    _validate_production_container_source()
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
    environment["PATH"] = os.pathsep.join(
        (str(Path(sys.executable).resolve().parent), environment.get("PATH", ""))
    )
    for variable in ("XDG_CACHE_HOME", "APPTAINER_CACHEDIR", "APPTAINER_TMPDIR"):
        Path(environment[variable]).mkdir(parents=True, exist_ok=True)

    print(
        "=== Fixture generation (purpose-built DB; uniprotDBRetrieve is not run) ===",
        flush=True,
    )
    _write_genome(workdir)
    fixture_snakefile = _write_fixture_snakefile(workdir)
    fixture_command = [
        snakemake,
        "fixture-build/prodigal.gff",
        "--snakefile",
        str(fixture_snakefile),
        "--directory",
        str(workdir),
        *_deployment_arguments(apptainer_prefix),
        "--software-deployment-method",
        "apptainer",
        "--allowed-rules",
        "fixtureProdigal",
    ]
    fixture_output = _stream(
        fixture_command,
        environment,
        workdir / "logs" / "harness-fixture-prodigal.log",
    )
    if "rule fixtureProdigal:" not in fixture_output:
        raise AssertionError("fixture Prodigal rule did not execute")
    _validate_cached_image(apptainer, apptainer_prefix)
    _finish_fixture(workdir, workdir / "fixture-build" / "prodigal.gff")
    input_snapshot = _input_snapshot(workdir)

    print(
        "=== Production rules (annotation adapter -> Reparation -> GFF) ===",
        flush=True,
    )
    production_snakefile = _write_production_snakefile(workdir)
    target = "reparation/A-1.reparation.gff"
    production_command = [
        snakemake,
        target,
        "--snakefile",
        str(production_snakefile),
        "--directory",
        str(workdir),
        *_deployment_arguments(apptainer_prefix),
        "--software-deployment-method",
        "conda",
        "apptainer",
        "--conda-prefix",
        str(conda_prefix),
        "--resources",
        "reparation_instances=1",
        "--allowed-rules",
        "prepareReparationAnnotation",
        "reparation",
        "reparationGFF",
    ]
    first = _stream(
        production_command,
        environment,
        workdir / "logs" / "harness-production-first.log",
    )
    for rule in ("prepareReparationAnnotation", "reparation", "reparationGFF"):
        if f"rule {rule}:" not in first:
            raise AssertionError(f"production rule did not execute: {rule}")
    if "rule uniprotDBRetrieve:" in first:
        raise AssertionError("production smoke unexpectedly scheduled the live UniProt download")
    cached_runner = any(
        "source-cache" in line and "run_reparation.py" in line
        for line in first.splitlines()
    )
    if "run_reparation.py (cached)" not in first or not cached_runner:
        raise AssertionError("container runner was not resolved through Snakemake source-cache")

    _validate_outputs(workdir, input_snapshot)
    before = _snapshot(workdir)
    second = _stream(
        production_command,
        environment,
        workdir / "logs" / "harness-production-noop.log",
    )
    if "Nothing to be done" not in second:
        raise AssertionError("second production Snakemake run was not a no-op")
    if _snapshot(workdir) != before:
        raise AssertionError("a no-op rerun changed a published Reparation artifact")

    print("Reparation Snakemake + Apptainer smoke passed", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (
        AssertionError,
        OSError,
        RuntimeError,
        ValueError,
        csv.Error,
        subprocess.SubprocessError,
    ) as error:
        print(f"Reparation container smoke failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
