#!/usr/bin/env python3
"""Stage user inputs without changing what downstream rules see.

FASTQ files remain symlinks, but their targets are made absolute so both
relative project paths and already-absolute paths work. Text references are
copied byte-for-byte unless their gzip magic bytes show that they need to be
decompressed first.
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import tempfile
from pathlib import Path

GZIP_MAGIC = b"\x1f\x8b"


def link_input(source: Path, destination: Path) -> None:
    """Atomically point ``destination`` at an existing absolute ``source``."""
    source = Path(os.path.abspath(os.fspath(source)))
    if not source.is_file():
        raise FileNotFoundError(f"input file does not exist or is not a file: {source}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination = Path(os.path.abspath(os.fspath(destination)))
    if source == destination:
        raise ValueError("input and staged-link destination must be different paths")

    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    temporary.unlink()

    try:
        temporary.symlink_to(source)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def materialize_text(source: Path, destination: Path) -> None:
    """Atomically copy ``source`` to plain text, decompressing gzip if needed."""
    destination.parent.mkdir(parents=True, exist_ok=True)

    with source.open("rb") as probe:
        compressed = probe.read(2) == GZIP_MAGIC

    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
    )
    os.close(descriptor)
    temporary = Path(temporary_name)

    try:
        opener = gzip.open if compressed else open
        with opener(source, "rb") as source_handle, temporary.open("wb") as output_handle:
            shutil.copyfileobj(source_handle, output_handle)
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create an input symlink or materialize a plain-text reference."
    )
    parser.add_argument("operation", choices=("link", "text"))
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    if args.operation == "link":
        link_input(args.source, args.destination)
    else:
        materialize_text(args.source, args.destination)


if __name__ == "__main__":
    main()
