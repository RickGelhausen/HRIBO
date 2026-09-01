#!/usr/bin/env python3
"""Concatenate GFF3 records into one deterministic, atomically written file."""

import argparse
import os
import sys
import tempfile
from pathlib import Path


GFF3_HEADER = "##gff-version 3\n"


def _coordinate(value, path, line_number):
    try:
        return int(value)
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number}: GFF3 coordinates must be integers, found {value!r}"
        ) from error


def read_records(paths):
    """Return validated GFF3 records from all inputs.

    Directives and comments are intentionally not copied. The aggregate owns one
    canonical version directive, and retaining per-input directives would put
    headers in the middle of the output.
    """

    records = []
    for filename in paths:
        path = Path(filename)
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.rstrip("\r\n")
                if not stripped or stripped.startswith("#"):
                    continue

                fields = stripped.split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"{path}:{line_number}: expected 9 tab-separated GFF3 columns, "
                        f"found {len(fields)}"
                    )

                start = _coordinate(fields[3], path, line_number)
                end = _coordinate(fields[4], path, line_number)
                records.append((fields, start, end))
    return records


def _feature_rank(feature):
    """Keep parents ahead of children when coordinates are otherwise identical."""

    feature = feature.lower()
    if feature in {"gene", "pseudogene"}:
        return 0
    if feature in {"mrna", "transcript"} or feature.endswith("rna"):
        return 1
    return 2


def sort_records(records):
    """Sort structurally, with complete deterministic tie breakers."""

    return sorted(
        records,
        key=lambda record: (
            record[0][0],
            record[1],
            record[2],
            record[0][6],
            _feature_rank(record[0][2]),
            record[0][2].lower(),
            record[0][1],
            record[0][5],
            record[0][7],
            record[0][8],
        ),
    )


def render(records):
    return GFF3_HEADER + "".join("\t".join(fields) + "\n" for fields, _, _ in records)


def atomic_write(path, content):
    """Replace *path* only after its complete new contents are on disk."""

    output = Path(path)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp", text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
        os.chmod(temporary_name, 0o644)
        os.replace(temporary_name, output)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def concatenate(input_files, output_file=""):
    content = render(sort_records(read_records(input_files)))
    if output_file:
        atomic_write(output_file, content)
    else:
        sys.stdout.write(content)


def main():
    parser = argparse.ArgumentParser(description="Concatenate GFF3 files.")
    parser.add_argument(
        "input_files", nargs="*", metavar="GFF3", help="GFF3 files to concatenate"
    )
    parser.add_argument(
        "-o",
        "--output-file",
        "--output_file",
        default="",
        help="output GFF3 (default: stdout)",
    )
    args = parser.parse_args()

    try:
        concatenate(args.input_files, args.output_file)
    except (OSError, ValueError) as error:
        parser.exit(1, f"concatenate_gff.py: error: {error}\n")


if __name__ == "__main__":
    main()
