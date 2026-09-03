#!/usr/bin/env python3
"""Run the pinned DeepRibo parser with a strict headerless-bedGraph loader.

The DeepRibo v1.1 parser calls ``read_csv(..., skiprows=1)`` even though both
its example files and HRIBO's generated bedGraphs have no header.  The parser
installed in HRIBO's pinned image additionally turns one-row contigs into
all-zero signal.  This launcher verifies that exact installed source, replaces
only ``loadSignal``, and then delegates every other operation to DeepRibo.

Added by the HRIBO project in 2026; distributed under GPL-3.0 with HRIBO and
DeepRibo.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import ModuleType

import numpy as np


IMAGE_PARSER = Path("/usr/local/bin/DataParser.py")
IMAGE_PARSER_SHA256 = (
    "b13827da386f4e88680b906df982edf744ebbbf5403c3e47bc257798c8e2590d"
)


class BedGraphError(ValueError):
    """A DeepRibo signal track is malformed or incompatible with the genome."""


@dataclass(frozen=True)
class BedGraphRow:
    chrom: str
    start: int
    end: int
    count: float
    line_number: int


@lru_cache(maxsize=8)
def read_bedgraph(path_string: str) -> tuple[BedGraphRow, ...]:
    """Read a headerless 0-based, half-open bedGraph, including empty files."""

    path = Path(path_string)
    rows = []
    previous_by_chrom: dict[str, BedGraphRow] = {}
    try:
        handle = path.open(encoding="utf-8")
    except OSError as exc:
        raise BedGraphError(f"cannot read {path}: {exc}") from exc

    with handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if (
                not stripped
                or stripped.startswith("#")
                or stripped.startswith("track ")
                or stripped.startswith("browser ")
            ):
                continue
            fields = stripped.split()
            if len(fields) != 4:
                raise BedGraphError(
                    f"{path}: line {line_number} has {len(fields)} fields; expected four"
                )
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
                count = float(fields[3])
            except ValueError as exc:
                raise BedGraphError(
                    f"{path}: line {line_number} has non-numeric coordinates or count"
                ) from exc
            if start < 0 or end <= start:
                raise BedGraphError(
                    f"{path}: line {line_number} has invalid half-open interval "
                    f"{start}-{end}"
                )
            if not math.isfinite(count) or count < 0:
                raise BedGraphError(
                    f"{path}: line {line_number} has invalid count {fields[3]!r}"
                )
            row = BedGraphRow(chrom, start, end, count, line_number)
            previous = previous_by_chrom.get(chrom)
            if previous is not None and start < previous.end:
                raise BedGraphError(
                    f"{path}: line {line_number} overlaps or precedes line "
                    f"{previous.line_number} on {chrom}"
                )
            previous_by_chrom[chrom] = row
            rows.append(row)
    return tuple(rows)


def load_signal(
    ribo_coverage: str,
    ribo_elongating: str,
    chrom: str,
    genome_length: int,
    asense: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Return complete signal vectors for one contig.

    ``asense`` is part of DeepRibo's public callback signature. Orientation is
    applied later by its parser, so the loader intentionally does not use it.
    """

    del asense

    def vector(path: str) -> np.ndarray:
        signal = np.zeros(genome_length)
        for row in read_bedgraph(str(Path(path).resolve())):
            if row.chrom != chrom:
                continue
            if row.end > genome_length:
                raise BedGraphError(
                    f"{path}: line {row.line_number} ends at {row.end}, beyond "
                    f"{chrom}'s length {genome_length}"
                )
            signal[row.start : row.end] = row.count
        return signal

    return vector(ribo_coverage), vector(ribo_elongating)


def fasta_lengths(path: Path) -> dict[str, int]:
    """Read the FASTA identifiers and lengths using DeepRibo's ID convention."""

    lengths: dict[str, int] = {}
    current = ""
    try:
        handle = path.open(encoding="utf-8")
    except OSError as exc:
        raise BedGraphError(f"cannot read {path}: {exc}") from exc

    with handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                identifier_fields = stripped[1:].split()
                if not identifier_fields:
                    raise BedGraphError(f"{path}: line {line_number} has no FASTA ID")
                current = identifier_fields[0]
                if current in lengths:
                    raise BedGraphError(f"{path}: duplicate FASTA ID {current!r}")
                lengths[current] = 0
            elif not current:
                raise BedGraphError(
                    f"{path}: sequence data precedes the first FASTA header"
                )
            else:
                lengths[current] += len(stripped)
    if not lengths:
        raise BedGraphError(f"{path} contains no FASTA records")
    return lengths


def validate_tracks(track_paths: list[str], fasta: str) -> None:
    """Reject unknown contigs and intervals beyond their FASTA sequences."""

    lengths = fasta_lengths(Path(fasta))
    for path_string in track_paths:
        for row in read_bedgraph(str(Path(path_string).resolve())):
            length = lengths.get(row.chrom)
            if length is None:
                raise BedGraphError(
                    f"{path_string}: line {row.line_number} names unknown contig "
                    f"{row.chrom!r}"
                )
            if row.end > length:
                raise BedGraphError(
                    f"{path_string}: line {row.line_number} ends at {row.end}, "
                    f"beyond {row.chrom}'s length {length}"
                )


def load_image_parser(path: Path = IMAGE_PARSER) -> ModuleType:
    """Load only the exact parser audited with HRIBO's image digest."""

    try:
        contents = path.read_bytes()
    except OSError as exc:
        raise RuntimeError(f"cannot read pinned DeepRibo parser {path}: {exc}") from exc
    actual = hashlib.sha256(contents).hexdigest()
    if actual != IMAGE_PARSER_SHA256:
        raise RuntimeError(
            f"unexpected DeepRibo parser at {path}: SHA-256 {actual}, "
            f"expected {IMAGE_PARSER_SHA256}"
        )
    spec = importlib.util.spec_from_file_location("hribo_pinned_deepribo", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot create an import specification for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    try:
        # Leave usage handling to DeepRibo's own argparse parser. When all five
        # positional input paths are present, validate them before it writes any
        # tensors or data_list.csv rows.
        if len(sys.argv) >= 6:
            validate_tracks(sys.argv[1:5], sys.argv[5])
        parser = load_image_parser()
        parser.loadSignal = load_signal
        parser.main()
    except (BedGraphError, RuntimeError) as exc:
        raise SystemExit(f"deepribo_data_parser: error: {exc}") from exc


if __name__ == "__main__":
    main()
