#!/usr/bin/env python3
"""Compare primary outputs from two HRIBO real-data runs.

The comparator is intentionally semantic.  A release-validation run should not
require byte-identical scientific output after bug fixes or dependency updates,
but it should make missing artifacts, malformed files, changed feature sets and
large numerical shifts visible in one deterministic report.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import json
import math
import os
import statistics
import struct
import sys
import tempfile
from collections import Counter, defaultdict, deque
from pathlib import Path

import openpyxl
import pysam


REPORT_SCHEMA_VERSION = 1

# These are stable or scientifically useful outputs shared by HRIBO 1.8 and the
# 2.0 workflow.  Plots are deliberately not compared pixel-by-pixel; their
# source tables are included instead.
OUTPUT_PATTERNS = (
    "auxiliary/*.xlsx",
    "auxiliary/overview.tsv",
    "auxiliary/overview.gff",
    "auxiliary/overview_misc.gff",
    "metageneprofiling/**/*.xlsx",
    "readcounts/*.csv",
    "pca/*.csv",
    "pca/*.tsv",
    "xtail/*.csv",
    "xtail/*_sorted.xlsx",
    "riborex/*.csv",
    "riborex/*_sorted.xlsx",
    "deltate/*.csv",
    "deltate/*_sorted.xlsx",
    "figures/SpearmanCorr_readCounts.tab",
    "SpearmanCorr_readCounts.tab",
    "globaltracks/**/*.bw",
    "globaltracks/**/*.bigWig",
    "globaltracks/**/*.bigwig",
    "centeredtracks/**/*.bw",
    "centeredtracks/**/*.bigWig",
    "centeredtracks/**/*.bigwig",
    "fiveprimetracks/**/*.bw",
    "fiveprimetracks/**/*.bigWig",
    "fiveprimetracks/**/*.bigwig",
    "threeprimetracks/**/*.bw",
    "threeprimetracks/**/*.bigWig",
    "threeprimetracks/**/*.bigwig",
    "tracks/*.gff",
    "maplink/*.bam",
    "tis_advice/**/*.json",
    "tis_advice/**/*.tsv",
)

KEY_COLUMN_CANDIDATES = (
    ("identifier", "contrast"),
    ("method", "condition", "replicate"),
    ("identifier", "feature"),
    ("identifier",),
    ("genome", "start", "stop", "strand"),
    ("seqid", "start", "stop", "strand"),
    ("orientation", "class"),
    ("position",),
    ("reference",),
    ("sample",),
    ("id",),
)

HEADER_ALIASES = {
    "identifer": "identifier",
    "gene_id": "identifier",
    "log2foldchange": "log2fc",
    "lfcse": "log2fc_se",
    "log2foldchange_se": "log2fc_se",
    "padj": "pvalue_adjusted",
    "pvalue_adjust": "pvalue_adjusted",
    "adjusted_pvalue": "pvalue_adjusted",
    "pred_value": "deepribo_score",
    "pred_rank": "deepribo_rank",
    "pred_probability": "reparation_probability",
    "15nt_upstream": "upstream_15nt",
}

MISSING_TEXT_TOKENS = {
    "<na>",
    "#n/a",
    "n/a",
    "n.a.",
    "na",
    "na_character_",
    "na_complex_",
    "na_integer_",
    "na_real_",
    "none",
    "null",
}

BIGWIG_MAGIC = 0x888FFC26
BIGWIG_CHROM_TREE_MAGIC = 0x78CA8C91
LEGACY_DEEPRIBO_GFF_NAMES = frozenset(
    {
        "deepribo_all.gff",
        "deepribo_merged.gff",
        "deepribo_merged_plus.gff",
        "totalannotation.gff",
        "updated_annotation.gff",
    }
)
LEGACY_EMPTY_CONDITION_GFF_SUFFIXES = (
    ".deepribo.gff",
    ".reparation.gff",
    ".merged.gff",
)


class ComparisonError(ValueError):
    """Raised when a run or artifact cannot be inspected safely."""


def resolve_path(path: Path, *, description: str, strict: bool = False) -> Path:
    """Resolve a path while translating symlink failures into a clean error."""

    try:
        return path.resolve(strict=strict)
    except (OSError, RuntimeError) as error:
        raise ComparisonError(f"cannot resolve {description}: {error}") from error


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest for a regular-sized result file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_header(value: object) -> str:
    """Normalize harmless historical header spelling differences."""

    if value is None:
        return ""
    text = str(value).strip().lower()
    result = []
    previous_separator = False
    for character in text:
        if character.isalnum():
            result.append(character)
            previous_separator = False
        elif not previous_separator:
            result.append("_")
            previous_separator = True
    normalized = "".join(result).strip("_")
    if normalized in HEADER_ALIASES:
        return HEADER_ALIASES[normalized]
    # deltaTE prefixes the same DESeq2 column names with RIBO/RNA/TE.  Apply
    # aliases to that suffix too so old and current workbooks line up.
    for legacy, current in (
        ("log2foldchange", "log2fc"),
        ("log2foldchange_se", "log2fc_se"),
        ("lfcse", "log2fc_se"),
        ("padj", "pvalue_adjusted"),
        ("pvalue_adjust", "pvalue_adjusted"),
    ):
        suffix = f"_{legacy}"
        if normalized.endswith(suffix):
            return f"{normalized[:-len(suffix)]}_{current}"
    return normalized


def display_value(value: object) -> str:
    """Create a stable, compact value representation for keys and examples."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return format(value, ".15g")
    return str(value).strip()


def numeric_value(value: object) -> float | None:
    """Return a finite float for numeric cells, otherwise ``None``."""

    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        number = float(value)
    elif isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        try:
            number = float(stripped)
        except ValueError:
            return None
    else:
        return None
    return number if math.isfinite(number) else None


def nonfinite_numeric_token(value: object) -> bool:
    """Recognize explicit NaN/infinity values without treating textual NA as one."""

    if isinstance(value, float):
        return not math.isfinite(value)
    if not isinstance(value, str):
        return False
    return value.strip().lower() in {
        "nan",
        "+nan",
        "-nan",
        "inf",
        "+inf",
        "-inf",
        "infinity",
        "+infinity",
        "-infinity",
    }


def missing_value(value: object) -> bool:
    if value is None:
        return True
    return isinstance(value, str) and (
        not value.strip() or value.strip().lower() in MISSING_TEXT_TOKENS
    )


def infer_key_columns(
    columns: tuple[str, ...],
    *,
    relative_path: str = "",
) -> tuple[str, ...]:
    """Choose the most stable available key, falling back to the first column."""

    available = set(columns)
    lower_path = relative_path.lower()
    if lower_path == "pca/rld.tsv" and {"sampletype", "name"} <= available:
        return ("sampletype", "name")
    if lower_path in {
        "pca/meta.csv",
        "pca/rld_cor.tsv",
        "figures/spearmancorr_readcounts.tab",
    } and "sample" in available:
        return ("sample",)
    for candidate in KEY_COLUMN_CANDIDATES:
        if set(candidate) <= available:
            return candidate
    return columns[:1]


def make_key(row: dict[str, object], columns: tuple[str, ...]) -> str:
    values = [display_value(row.get(column)) for column in columns]
    return json.dumps(values, ensure_ascii=True, separators=(",", ":"))


def build_table(
    headers: list[object],
    raw_rows: list[list[object]],
    *,
    relative_path: str = "",
    strict_row_shape: bool = False,
) -> dict:
    """Build an internal normalized table plus a JSON-safe public summary."""

    original_columns = tuple(display_value(header) for header in headers)
    columns = tuple(normalize_header(header) for header in headers)
    if not columns or any(not column for column in columns):
        raise ComparisonError("table has an empty column heading")
    duplicates = sorted(
        column for column, count in Counter(columns).items() if count > 1
    )
    if duplicates:
        raise ComparisonError(
            "table has duplicate normalized column headings: " + ", ".join(duplicates)
        )

    rows = []
    for row_number, raw_row in enumerate(raw_rows, start=2):
        if not any(display_value(value) for value in raw_row):
            continue
        if strict_row_shape and len(raw_row) != len(columns):
            raise ComparisonError(
                f"table row {row_number} has {len(raw_row)} columns instead of "
                f"{len(columns)}"
            )
        if len(raw_row) > len(columns) and any(
            display_value(value) for value in raw_row[len(columns) :]
        ):
            raise ComparisonError(
                f"table row {row_number} has populated cells beyond its "
                f"{len(columns)} headings"
            )
        padded = list(raw_row[: len(columns)])
        padded.extend([None] * (len(columns) - len(padded)))
        rows.append(dict(zip(columns, padded, strict=True)))

    key_columns = infer_key_columns(columns, relative_path=relative_path)
    keyed_rows: dict[str, dict[str, object]] = {}
    duplicate_keys = []
    for index, row in enumerate(rows):
        key = make_key(row, key_columns)
        if not any(display_value(row.get(column)) for column in key_columns):
            key = json.dumps(["__row__", index], separators=(",", ":"))
        if key in keyed_rows:
            duplicate_keys.append(key)
            # Keep every row comparable without silently overwriting a duplicate.
            key = json.dumps(["__duplicate__", key, index], separators=(",", ":"))
        keyed_rows[key] = row

    numeric_columns = []
    nonfinite_cells = 0
    missing_cells = 0
    for column in columns:
        missing_cells += sum(missing_value(row[column]) for row in rows)
        populated = [row[column] for row in rows if not missing_value(row[column])]
        nonfinite_cells += sum(nonfinite_numeric_token(value) for value in populated)
        if populated and all(
            numeric_value(value) is not None or nonfinite_numeric_token(value)
            for value in populated
        ):
            numeric_columns.append(column)

    return {
        "columns": columns,
        "original_columns": original_columns,
        "rows": rows,
        "key_columns": key_columns,
        "keyed_rows": keyed_rows,
        "summary": {
            "rows": len(rows),
            "columns": list(original_columns),
            "normalized_columns": list(columns),
            "key_columns": list(key_columns),
            "duplicate_keys": len(duplicate_keys),
            "duplicate_key_examples": duplicate_keys[:5],
            "numeric_columns": numeric_columns,
            "missing_cells": missing_cells,
            "nonfinite_numeric_cells": nonfinite_cells,
        },
    }


def read_delimited_table(path: Path, relative_path: str) -> dict[str, dict]:
    """Read one of HRIBO's delimited formats using its actual on-disk dialect."""

    lower_path = relative_path.lower()
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    row_heading = None
    headerless = False
    if lower_path == "pca/meta.csv":
        delimiter = "\t"  # Historical extension; preparePCAinput writes TSV.
        row_heading = "Sample"
    elif lower_path == "pca/rld_cor.tsv":
        row_heading = "Sample"
    elif lower_path == "pca/variance_percentages.tsv":
        headerless = True
    elif lower_path == "figures/spearmancorr_readcounts.tab":
        row_heading = "Sample"
    elif lower_path.startswith(("xtail/", "riborex/")):
        # R's write.csv defaults to row.names=TRUE and leaves this heading blank.
        row_heading = "Identifier"

    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = [
                list(row)
                for row in csv.reader(handle, delimiter=delimiter, strict=True)
            ]
    except (OSError, UnicodeError, csv.Error) as error:
        raise ComparisonError(f"cannot parse delimited table: {error}") from error

    rows = [row for row in rows if any(display_value(value) for value in row)]
    if lower_path == "figures/spearmancorr_readcounts.tab":
        while rows and len(rows[0]) == 1 and rows[0][0].startswith("#plotCorrelation"):
            rows.pop(0)
    if not rows:
        raise ComparisonError("table is empty")

    if headerless:
        for line_number, row in enumerate(rows, start=1):
            if len(row) != 1:
                raise ComparisonError(
                    f"table row {line_number} has {len(row)} columns instead of 1"
                )
        data = [[f"PC{index}", row[0]] for index, row in enumerate(rows, start=1)]
        return {
            "table": build_table(
                ["Component", "Variance percentage"],
                data,
                relative_path=relative_path,
                strict_row_shape=True,
            )
        }

    headers, data = rows[0], rows[1:]
    if row_heading is not None:
        if headers and not display_value(headers[0]):
            headers[0] = row_heading
        elif data and all(len(row) == len(headers) + 1 for row in data):
            # Some R/deepTools matrix writers omit the row-name placeholder.
            headers.insert(0, row_heading)
    return {
        "table": build_table(
            headers,
            data,
            relative_path=relative_path,
            strict_row_shape=True,
        )
    }


def read_workbook(path: Path, relative_path: str) -> dict[str, dict]:
    try:
        workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    except Exception as error:
        raise ComparisonError(f"cannot open workbook: {error}") from error

    tables = {}
    try:
        for worksheet in workbook.worksheets:
            row_iterator = worksheet.iter_rows(values_only=True)
            headers = None
            raw_rows = []
            for values in row_iterator:
                values = list(values)
                if headers is None:
                    if any(display_value(value) for value in values):
                        headers = values
                    continue
                raw_rows.append(values)
            if headers is None:
                raise ComparisonError(f"worksheet {worksheet.title!r} is empty")
            tables[worksheet.title] = build_table(
                headers, raw_rows, relative_path=relative_path
            )
    finally:
        workbook.close()
    if not tables:
        raise ComparisonError("workbook has no worksheets")
    return tables


def parse_gff_attributes(text: str) -> dict[str, str]:
    text = text.strip()
    if text == ".":
        return {}
    if not text:
        raise ComparisonError("GFF attributes column is empty; use '.' when absent")
    attributes = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ComparisonError(f"GFF attribute token {item!r} is not key=value")
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise ComparisonError("GFF attribute has an empty key")
        if key in attributes:
            raise ComparisonError(f"GFF attribute key {key!r} occurs more than once")
        attributes[key] = value
    if not attributes:
        raise ComparisonError("GFF attributes contain no key=value token; use '.'")
    return attributes


def is_legacy_deepribo_gff_path(relative_path: str) -> bool:
    """Return whether HRIBO 1.8 could publish DeepRibo phase misuse here."""

    if not relative_path.startswith("tracks/") or relative_path.count("/") != 1:
        return False
    filename = relative_path.rsplit("/", 1)[1].lower()
    return filename in LEGACY_DEEPRIBO_GFF_NAMES or filename.endswith(
        ".deepribo.gff"
    )


def is_legacy_empty_condition_gff_path(relative_path: str) -> bool:
    """Return whether HRIBO 1.8 could emit an empty per-condition GFF here."""

    if not relative_path.startswith("tracks/") or relative_path.count("/") != 1:
        return False
    filename = relative_path.rsplit("/", 1)[1].lower()
    return any(
        filename.endswith(suffix) and len(filename) > len(suffix)
        for suffix in LEGACY_EMPTY_CONDITION_GFF_SUFFIXES
    )


def read_gff(
    path: Path,
    *,
    allow_legacy_deepribo_phase: bool = False,
    allow_empty_without_version: bool = False,
) -> dict:
    records = []
    legacy_deepribo_phases = []
    declares_gff3 = False
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if line == "##gff-version 3":
                    declares_gff3 = True
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) != 9:
                    raise ComparisonError(
                        f"GFF line {line_number} has {len(fields)} columns instead of 9"
                    )
                if not fields[0].strip() or fields[0].strip() == ".":
                    raise ComparisonError(
                        f"GFF line {line_number} has an empty/placeholder seqid"
                    )
                if not fields[1].strip():
                    raise ComparisonError(f"GFF line {line_number} has an empty source")
                if not fields[2].strip() or fields[2].strip() == ".":
                    raise ComparisonError(
                        f"GFF line {line_number} has an empty/placeholder feature"
                    )
                try:
                    start = int(fields[3])
                    end = int(fields[4])
                except ValueError as error:
                    raise ComparisonError(
                        f"GFF line {line_number} has non-integer coordinates"
                    ) from error
                if start < 1 or end < start:
                    raise ComparisonError(
                        f"GFF line {line_number} has invalid coordinates {start}-{end}"
                    )
                if fields[6] not in {"+", "-", ".", "?"}:
                    raise ComparisonError(
                        f"GFF line {line_number} has invalid strand {fields[6]!r}"
                    )
                try:
                    attributes = parse_gff_attributes(fields[8])
                except ComparisonError as error:
                    raise ComparisonError(f"GFF line {line_number}: {error}") from error
                score_value = None
                if fields[5] != ".":
                    try:
                        score_value = float(fields[5])
                    except ValueError as error:
                        raise ComparisonError(
                            f"GFF line {line_number} has a non-numeric score"
                        ) from error
                    if not math.isfinite(score_value):
                        raise ComparisonError(
                            f"GFF line {line_number} has a non-finite score"
                        )
                phase = fields[7]
                attribute_keys = {key.lower() for key in attributes}
                legacy_deepribo_phase = bool(
                    allow_legacy_deepribo_phase
                    and fields[1].strip().lower() == "deepribo"
                    and fields[2].strip().lower() == "cds"
                    and "deepribo_distance" not in attribute_keys
                    and phase != "."
                    and phase.removeprefix("-").isascii()
                    and phase.removeprefix("-").isdigit()
                )
                if legacy_deepribo_phase:
                    # HRIBO 1.8 stored DeepRibo distance/novel-rank metadata in
                    # the GFF3 phase field.  Current output moves that metadata
                    # to attributes and uses the required zero CDS phase.
                    legacy_deepribo_phases.append(phase)
                    phase = "0"
                elif phase not in {".", "0", "1", "2"}:
                    raise ComparisonError(
                        f"GFF line {line_number} has invalid phase {phase!r}"
                    )
                records.append(
                    {
                        "seqid": fields[0],
                        "source": fields[1],
                        "feature": fields[2],
                        "start": start,
                        "end": end,
                        "strand": fields[6],
                        "score": score_value,
                        "phase": phase,
                        "attributes": attributes,
                    }
                )
    except (OSError, UnicodeError) as error:
        raise ComparisonError(f"cannot parse GFF: {error}") from error

    if not records and not declares_gff3 and not allow_empty_without_version:
        raise ComparisonError(
            "GFF contains no feature records and does not declare "
            "'##gff-version 3'"
        )

    keys = [gff_key(record) for record in records]
    summary = {
        "records": len(records),
        "duplicate_coordinate_keys": len(keys) - len(set(keys)),
        "by_feature": dict(sorted(Counter(r["feature"] for r in records).items())),
        "by_source": dict(sorted(Counter(r["source"] for r in records).items())),
        "by_strand": dict(sorted(Counter(r["strand"] for r in records).items())),
        "contigs": sorted({r["seqid"] for r in records}),
    }
    if legacy_deepribo_phases:
        summary["legacy_deepribo_phase_compatibility"] = {
            "records_normalized": len(legacy_deepribo_phases),
            "original_phase_counts": dict(
                sorted(Counter(legacy_deepribo_phases).items())
            ),
            "normalized_phase": "0",
        }
    return {
        "records": records,
        "keys": keys,
        "summary": summary,
    }


def gff_key(record: dict) -> tuple:
    return (
        record["seqid"],
        record["feature"],
        record["start"],
        record["end"],
        record["strand"],
    )


def gff_semantic_key(record: dict) -> tuple:
    """Represent a GFF record without depending on line/attribute order."""

    return (
        record["seqid"],
        record["source"],
        record["feature"],
        record["start"],
        record["end"],
        record["score"],
        record["strand"],
        record["phase"],
        tuple(sorted(record["attributes"].items())),
    )


def read_binary_region(
    handle, *, offset: int, length: int, file_size: int, label: str
) -> bytes:
    if offset < 0 or length < 0 or offset + length > file_size:
        raise ComparisonError(f"BigWig {label} points outside the file")
    handle.seek(offset)
    data = handle.read(length)
    if len(data) != length:
        raise ComparisonError(f"BigWig {label} is truncated")
    return data


def read_bigwig(path: Path) -> dict:
    """Validate and summarize the dependency-free portion of a BigWig file.

    Interval blocks can be compressed and indexed in several ways, so without
    pyBigWig the conservative equality contract remains the compressed-file
    digest.  The fixed header, zoom metadata, chromosome B+ tree and global
    numeric summary are still decoded and exposed for meaningful review.
    """

    file_size = path.stat().st_size
    try:
        with path.open("rb") as handle:
            magic_bytes = read_binary_region(
                handle, offset=0, length=4, file_size=file_size, label="magic"
            )
            if struct.unpack("<I", magic_bytes)[0] == BIGWIG_MAGIC:
                endian = "<"
                byte_order = "little"
            elif struct.unpack(">I", magic_bytes)[0] == BIGWIG_MAGIC:
                endian = ">"
                byte_order = "big"
            else:
                raise ComparisonError("file does not have the BigWig magic number")

            header = struct.unpack(
                f"{endian}IHHQQQHHQQIQ",
                read_binary_region(
                    handle,
                    offset=0,
                    length=64,
                    file_size=file_size,
                    label="fixed header",
                ),
            )
            (
                _magic,
                version,
                zoom_count,
                chromosome_tree_offset,
                full_data_offset,
                full_index_offset,
                field_count,
                defined_field_count,
                auto_sql_offset,
                total_summary_offset,
                uncompress_buffer_size,
                _reserved,
            ) = header
            if not 1 <= version <= 4:
                raise ComparisonError(f"BigWig has unsupported version {version}")
            if 64 + zoom_count * 24 > file_size:
                raise ComparisonError("BigWig zoom headers are truncated")
            for name, offset in (
                ("chromosome tree", chromosome_tree_offset),
                ("full data", full_data_offset),
                ("full index", full_index_offset),
            ):
                if offset < 64 or offset >= file_size:
                    raise ComparisonError(f"BigWig {name} offset is outside the file")
            for name, offset in (
                ("autoSql", auto_sql_offset),
                ("total summary", total_summary_offset),
            ):
                if offset and not 64 <= offset < file_size:
                    raise ComparisonError(f"BigWig {name} offset is outside the file")

            zoom_levels = []
            for index in range(zoom_count):
                reduction, _reserved, data_offset, index_offset = struct.unpack(
                    f"{endian}IIQQ",
                    read_binary_region(
                        handle,
                        offset=64 + index * 24,
                        length=24,
                        file_size=file_size,
                        label=f"zoom header {index}",
                    ),
                )
                for name, offset in (("data", data_offset), ("index", index_offset)):
                    if offset < 64 or offset >= file_size:
                        raise ComparisonError(
                            f"BigWig zoom {index} {name} offset is outside the file"
                        )
                zoom_levels.append(
                    {
                        "reduction": reduction,
                        "data_offset": data_offset,
                        "index_offset": index_offset,
                    }
                )

            tree_header = struct.unpack(
                f"{endian}IIIIQQ",
                read_binary_region(
                    handle,
                    offset=chromosome_tree_offset,
                    length=32,
                    file_size=file_size,
                    label="chromosome-tree header",
                ),
            )
            tree_magic, block_size, key_size, value_size, item_count, _reserved = (
                tree_header
            )
            if tree_magic != BIGWIG_CHROM_TREE_MAGIC:
                raise ComparisonError("BigWig chromosome tree has an invalid magic number")
            if not block_size or not key_size or value_size != 8:
                raise ComparisonError("BigWig chromosome tree has invalid dimensions")
            minimum_leaf_bytes = key_size + value_size
            if item_count > file_size // minimum_leaf_bytes:
                raise ComparisonError("BigWig chromosome tree item count is implausible")

            chromosomes = {}
            chromosome_ids = set()
            root_offset = chromosome_tree_offset + 32
            pending = [root_offset]
            visited = set()
            max_nodes = file_size // 4 + 1
            while pending:
                node_offset = pending.pop()
                if node_offset in visited:
                    raise ComparisonError("BigWig chromosome tree contains a cycle")
                visited.add(node_offset)
                if len(visited) > max_nodes:
                    raise ComparisonError("BigWig chromosome tree has too many nodes")
                leaf, reserved, count = struct.unpack(
                    f"{endian}BBH",
                    read_binary_region(
                        handle,
                        offset=node_offset,
                        length=4,
                        file_size=file_size,
                        label="chromosome-tree node",
                    ),
                )
                if leaf not in {0, 1} or reserved != 0 or count > block_size:
                    raise ComparisonError("BigWig chromosome tree node is malformed")
                entry_size = key_size + (8 if leaf else 8)
                node_data = read_binary_region(
                    handle,
                    offset=node_offset + 4,
                    length=count * entry_size,
                    file_size=file_size,
                    label="chromosome-tree entries",
                )
                for index in range(count):
                    entry = node_data[index * entry_size : (index + 1) * entry_size]
                    key_bytes = entry[:key_size].split(b"\0", 1)[0]
                    try:
                        key = key_bytes.decode("utf-8")
                    except UnicodeError as error:
                        raise ComparisonError(
                            "BigWig chromosome tree has a non-UTF-8 name"
                        ) from error
                    if not key:
                        raise ComparisonError("BigWig chromosome tree has an empty name")
                    payload = entry[key_size:]
                    if leaf:
                        chromosome_id, chromosome_size = struct.unpack(
                            f"{endian}II", payload
                        )
                        if not chromosome_size:
                            raise ComparisonError(
                                f"BigWig chromosome {key!r} has zero length"
                            )
                        if key in chromosomes or chromosome_id in chromosome_ids:
                            raise ComparisonError(
                                "BigWig chromosome tree has duplicate names or IDs"
                            )
                        chromosomes[key] = {
                            "id": chromosome_id,
                            "length": chromosome_size,
                        }
                        chromosome_ids.add(chromosome_id)
                    else:
                        (child_offset,) = struct.unpack(f"{endian}Q", payload)
                        if child_offset < root_offset or child_offset >= file_size:
                            raise ComparisonError(
                                "BigWig chromosome-tree child points outside the tree"
                            )
                        pending.append(child_offset)
            if len(chromosomes) != item_count:
                raise ComparisonError(
                    "BigWig chromosome-tree item count does not match its leaves"
                )

            total_summary = None
            if total_summary_offset:
                bases, minimum, maximum, total, squares = struct.unpack(
                    f"{endian}Qdddd",
                    read_binary_region(
                        handle,
                        offset=total_summary_offset,
                        length=40,
                        file_size=file_size,
                        label="total summary",
                    ),
                )
                values = (minimum, maximum, total, squares)
                if not all(math.isfinite(value) for value in values):
                    raise ComparisonError("BigWig total summary contains non-finite values")
                if bases and minimum > maximum:
                    raise ComparisonError("BigWig total summary minimum exceeds maximum")
                total_summary = {
                    "bases_covered": bases,
                    "minimum": minimum,
                    "maximum": maximum,
                    "sum": total,
                    "sum_squares": squares,
                }
    except ComparisonError:
        raise
    except (OSError, struct.error) as error:
        raise ComparisonError(f"cannot parse BigWig: {error}") from error

    return {
        "summary": {
            "comparison_mode": (
                "validated header/chromosomes/global summary; conservative "
                "compressed-byte equality"
            ),
            "byte_order": byte_order,
            "version": version,
            "zoom_levels": zoom_levels,
            "field_count": field_count,
            "defined_field_count": defined_field_count,
            "uncompress_buffer_size": uncompress_buffer_size,
            "chromosomes": dict(sorted(chromosomes.items())),
            "total_summary": total_summary,
        }
    }


def tolerant_gff_matches(
    baseline: list[dict], candidate: list[dict], tolerance: int = 3
) -> int:
    """Return a deterministic maximum-cardinality one-to-one tolerance match."""

    if not baseline or not candidate:
        return 0
    # The common case avoids constructing a graph and preserves duplicates.
    if Counter(map(gff_key, baseline)) == Counter(map(gff_key, candidate)):
        return len(baseline)

    grouped_old: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    grouped_new: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for record in baseline:
        grouped_old[(record["seqid"], record["feature"], record["strand"])].append(
            record
        )
    for record in candidate:
        grouped_new[(record["seqid"], record["feature"], record["strand"])].append(
            record
        )

    matches = 0
    for group in sorted(set(grouped_old) & set(grouped_new)):
        old_records = sorted(grouped_old[group], key=gff_key)
        new_records = sorted(grouped_new[group], key=gff_key)
        by_start: dict[int, list[int]] = defaultdict(list)
        for new_index, record in enumerate(new_records):
            by_start[record["start"]].append(new_index)

        adjacency: list[list[int]] = []
        for old in old_records:
            options = []
            for start in range(old["start"] - tolerance, old["start"] + tolerance + 1):
                for new_index in by_start.get(start, ()):
                    new = new_records[new_index]
                    end_delta = abs(old["end"] - new["end"])
                    if end_delta <= tolerance:
                        start_delta = abs(old["start"] - new["start"])
                        options.append(
                            (
                                start_delta + end_delta,
                                start_delta,
                                end_delta,
                                gff_key(new),
                                new_index,
                            )
                        )
            adjacency.append([option[-1] for option in sorted(options)])
        matches += maximum_bipartite_matches(adjacency, len(new_records))
    return matches


def maximum_bipartite_matches(adjacency: list[list[int]], right_size: int) -> int:
    """Hopcroft-Karp cardinality for a deterministic, pre-sorted adjacency list."""

    left_size = len(adjacency)
    left_match = [-1] * left_size
    right_match = [-1] * right_size
    distance = [0] * left_size
    infinity = left_size + 1

    def breadth_first() -> bool:
        queue = deque()
        found = False
        for left in range(left_size):
            if left_match[left] == -1:
                distance[left] = 0
                queue.append(left)
            else:
                distance[left] = infinity
        while queue:
            left = queue.popleft()
            for right in adjacency[left]:
                paired = right_match[right]
                if paired == -1:
                    found = True
                elif distance[paired] == infinity:
                    distance[paired] = distance[left] + 1
                    queue.append(paired)
        return found

    def depth_first(start: int) -> bool:
        # Keep the Hopcroft-Karp search iterative: a long chain of nearby GFF
        # records must not hit Python's recursion limit.
        left_path = [start]
        next_option = [0]
        right_path: list[int] = []
        while left_path:
            left = left_path[-1]
            while next_option[-1] < len(adjacency[left]):
                right = adjacency[left][next_option[-1]]
                next_option[-1] += 1
                paired = right_match[right]
                if paired == -1:
                    for path_left, path_right in zip(
                        left_path, [*right_path, right], strict=True
                    ):
                        left_match[path_left] = path_right
                        right_match[path_right] = path_left
                    return True
                if distance[paired] == distance[left] + 1:
                    right_path.append(right)
                    left_path.append(paired)
                    next_option.append(0)
                    break
            else:
                distance[left] = infinity
                left_path.pop()
                next_option.pop()
                if right_path:
                    right_path.pop()
        return False

    count = 0
    while breadth_first():
        for left in range(left_size):
            if left_match[left] == -1 and depth_first(left):
                count += 1
    return count


def json_safe_value(value: object) -> object:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, bytes):
        return {"bytes_hex": value.hex()}
    if isinstance(value, (list, tuple)):
        return [json_safe_value(item) for item in value]
    if hasattr(value, "tolist"):
        return json_safe_value(value.tolist())
    return repr(value)


def alignment_multiset_digest(alignment: pysam.AlignmentFile) -> dict:
    """Hash every alignment while making the aggregate independent of BAM order."""

    modulus = 1 << 256
    digest_sum = 0
    digest_xor = 0
    count = 0
    for record in alignment.fetch(until_eof=True):
        tags = [
            (tag, value_type, json_safe_value(value))
            for tag, value, value_type in record.get_tags(with_value_type=True)
        ]
        tags.sort(
            key=lambda item: (
                item[0],
                item[1],
                json.dumps(item[2], sort_keys=True, ensure_ascii=True),
            )
        )
        fields = [
            record.query_name,
            record.flag,
            record.reference_name,
            record.reference_start,
            record.mapping_quality,
            record.cigarstring,
            record.next_reference_name,
            record.next_reference_start,
            record.template_length,
            record.query_sequence,
            list(record.query_qualities) if record.query_qualities is not None else None,
            tags,
        ]
        encoded = json.dumps(
            fields, ensure_ascii=True, separators=(",", ":"), sort_keys=True
        ).encode("utf-8")
        integer = int.from_bytes(hashlib.sha256(encoded).digest(), "big")
        digest_sum = (digest_sum + integer) % modulus
        digest_xor ^= integer
        count += 1
    return {
        "algorithm": "sha256-multiset-sum-xor-v1",
        "alignments": count,
        "sum": f"{digest_sum:064x}",
        "xor": f"{digest_xor:064x}",
    }


def read_bam(path: Path) -> dict:
    try:
        with pysam.AlignmentFile(str(path), "rb") as alignment:
            if not alignment.has_index():
                raise ComparisonError("BAM has no readable index")
            references = list(alignment.references)
            lengths = list(alignment.lengths)
            header = alignment.header.to_dict()
            index_stats = {item.contig: item for item in alignment.get_index_statistics()}
            rows = []
            for reference, length in zip(references, lengths, strict=True):
                stats = index_stats.get(reference)
                rows.append(
                    [
                        reference,
                        length,
                        stats.mapped if stats else 0,
                        stats.unmapped if stats else 0,
                    ]
                )
            alignment_digest = alignment_multiset_digest(alignment)
    except ComparisonError:
        raise
    except Exception as error:
        raise ComparisonError(f"cannot open indexed BAM: {error}") from error
    table = build_table(["Reference", "Length", "Mapped", "Unmapped"], rows)
    stable_header = {key: value for key, value in header.items() if key != "PG"}
    provenance = header.get("PG", [])
    stable_header_digest = hashlib.sha256(
        json.dumps(stable_header, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    ).hexdigest()
    provenance_digest = hashlib.sha256(
        json.dumps(provenance, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return {
        "tables": {"index_statistics": table},
        "semantic_signature": {
            "stable_header_sha256": stable_header_digest,
            "alignment_digest": alignment_digest,
        },
        "provenance_sha256": provenance_digest,
        "summary": {
            "comparison_mode": "order-independent full-alignment digest",
            "stable_header_sha256": stable_header_digest,
            "program_provenance_sha256": provenance_digest,
            "alignment_digest": alignment_digest,
            "references": len(references),
            "mapped": sum(row[2] for row in rows),
            "unmapped": sum(row[3] for row in rows),
            "per_reference": {
                row[0]: {
                    "length": row[1],
                    "mapped": row[2],
                    "unmapped": row[3],
                }
                for row in rows
            },
        },
    }


def reject_json_constant(token: str) -> None:
    raise ValueError(f"non-finite number {token!r}")


def read_json(path: Path) -> dict:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle, parse_constant=reject_json_constant)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ComparisonError(f"cannot parse JSON: {error}") from error
    return {
        "value": value,
        "summary": {
            "top_level_type": type(value).__name__,
            "top_level_keys": sorted(value) if isinstance(value, dict) else [],
        },
    }


def inspect_artifact(
    path: Path, relative_path: str = "", *, baseline_compatibility: bool = False
) -> dict:
    if not path.is_file():
        raise ComparisonError("path is missing, dangling, or not a regular file")
    file_size = path.stat().st_size
    suffix = path.suffix.lower()
    legacy_empty_condition_gff = bool(
        file_size == 0
        and suffix == ".gff"
        and baseline_compatibility
        and is_legacy_empty_condition_gff_path(relative_path)
    )
    if file_size == 0 and not legacy_empty_condition_gff:
        raise ComparisonError("file is empty")

    common = {"bytes": file_size}
    if suffix == ".xlsx":
        return {
            **common,
            "kind": "xlsx",
            "sha256": sha256_file(path),
            "tables": read_workbook(path, relative_path),
        }
    if suffix in {".csv", ".tsv", ".tab"}:
        return {
            **common,
            "kind": suffix[1:],
            "sha256": sha256_file(path),
            "tables": read_delimited_table(path, relative_path),
        }
    if suffix == ".gff":
        parsed = read_gff(
            path,
            allow_legacy_deepribo_phase=(
                baseline_compatibility
                and is_legacy_deepribo_gff_path(relative_path)
            ),
            allow_empty_without_version=legacy_empty_condition_gff,
        )
        if legacy_empty_condition_gff:
            parsed["summary"]["legacy_zero_byte_gff_compatibility"] = {
                "records_normalized": 0,
                "interpretation": "empty feature set",
            }
        return {
            **common,
            "kind": "gff",
            "sha256": sha256_file(path),
            **parsed,
        }
    if suffix == ".bam":
        return {
            **common,
            "kind": "bam",
            **read_bam(path),
        }
    if suffix in {".bw", ".bigwig"}:
        parsed = read_bigwig(path)
        return {
            **common,
            "kind": "bigwig",
            "sha256": sha256_file(path),
            **parsed,
        }
    if suffix == ".json":
        return {
            **common,
            "kind": "json",
            "sha256": sha256_file(path),
            **read_json(path),
        }
    raise ComparisonError(f"unsupported result type {suffix!r}")


def bam_index_sidecars(path: Path) -> tuple[Path, ...]:
    """Return every local BAM-index name recognized by htslib."""

    candidates = (
        Path(f"{path}.bai"),
        path.with_suffix(".bai"),
        Path(f"{path}.csi"),
        path.with_suffix(".csi"),
    )
    return tuple(dict.fromkeys(candidates))


def validate_bam_index_locations(path: Path, root: Path) -> None:
    """Reject BAM index sidecars whose symlink resolution escapes a run root."""

    for index_path in bam_index_sidecars(path):
        if not os.path.lexists(index_path):
            continue
        relative = index_path.relative_to(root).as_posix()
        resolved = resolve_path(
            index_path,
            description=f"BAM index {relative!r}",
            strict=True,
        )
        if not resolved.is_relative_to(root):
            raise ComparisonError(
                f"BAM index symlink escapes run root: {relative} -> {resolved}"
            )


def discover_outputs(root: Path) -> dict[str, Path]:
    root = resolve_path(root, description=f"run root {root}")
    if not root.is_dir():
        raise ComparisonError(f"run root is not a directory: {root}")
    outputs = {}
    for pattern in OUTPUT_PATTERNS:
        for path in root.glob(pattern):
            relative = path.relative_to(root).as_posix()
            resolved = resolve_path(
                path, description=f"output {relative!r}", strict=True
            )
            if not resolved.is_relative_to(root):
                raise ComparisonError(
                    f"output symlink escapes run root: {relative} -> {resolved}"
                )
            if path.suffix.lower() == ".bam":
                validate_bam_index_locations(path, root)
            canonical = (
                "figures/SpearmanCorr_readCounts.tab"
                if relative == "SpearmanCorr_readCounts.tab"
                else relative
            )
            if canonical in outputs and outputs[canonical] != path:
                raise ComparisonError(
                    f"multiple outputs map to canonical path {canonical!r}"
                )
            # Keep the validated logical path.  maplink BAMs are commonly
            # symlinks whose index exists beside the link, not beside its target.
            outputs[canonical] = path
    return dict(sorted(outputs.items()))


def average_ranks(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    cursor = 0
    while cursor < len(indexed):
        end = cursor + 1
        while end < len(indexed) and indexed[end][1] == indexed[cursor][1]:
            end += 1
        rank = (cursor + 1 + end) / 2
        for index, _value in indexed[cursor:end]:
            ranks[index] = rank
        cursor = end
    return ranks


def pearson(values_a: list[float], values_b: list[float]) -> float | None:
    if len(values_a) < 2 or len(values_a) != len(values_b):
        return None
    mean_a = statistics.fmean(values_a)
    mean_b = statistics.fmean(values_b)
    delta_a = [value - mean_a for value in values_a]
    delta_b = [value - mean_b for value in values_b]
    denominator = math.sqrt(
        sum(value * value for value in delta_a)
        * sum(value * value for value in delta_b)
    )
    if denominator == 0:
        return None
    return sum(a * b for a, b in zip(delta_a, delta_b, strict=True)) / denominator


def safe_round(value: float | None) -> float | None:
    return None if value is None or not math.isfinite(value) else round(value, 10)


def compare_numeric_column(
    old_rows: dict[str, dict[str, object]],
    new_rows: dict[str, dict[str, object]],
    keys: list[str],
    old_column: str,
    new_column: str,
) -> dict:
    pairs = []
    missing_both = 0
    nonfinite_pairs = 0
    nonfinite_value_mismatches = 0
    text_pairs = 0
    text_value_mismatches = 0
    missing_mismatches = 0
    nonfinite_mismatches = 0
    type_mismatches = 0
    mismatch_examples = []
    for key in keys:
        old_raw = old_rows[key].get(old_column)
        new_raw = new_rows[key].get(new_column)
        old_kind, old_value = classify_numeric_cell(old_raw)
        new_kind, new_value = classify_numeric_cell(new_raw)
        if old_kind == new_kind == "finite":
            pairs.append((old_value, new_value))
            continue
        if old_kind == new_kind == "missing":
            missing_both += 1
            continue
        mismatch = False
        if old_kind == new_kind == "nonfinite":
            nonfinite_pairs += 1
            mismatch = old_value != new_value
            nonfinite_value_mismatches += mismatch
        elif old_kind == new_kind == "text":
            text_pairs += 1
            mismatch = display_value(old_raw) != display_value(new_raw)
            text_value_mismatches += mismatch
        elif "missing" in {old_kind, new_kind}:
            missing_mismatches += 1
            mismatch = True
        elif "nonfinite" in {old_kind, new_kind}:
            nonfinite_mismatches += 1
            mismatch = True
        else:
            type_mismatches += 1
            mismatch = True
        if mismatch and len(mismatch_examples) < 10:
            mismatch_examples.append(
                {
                    "key": key,
                    "baseline": display_value(old_raw),
                    "candidate": display_value(new_raw),
                    "baseline_type": old_kind,
                    "candidate_type": new_kind,
                }
            )

    values_a = [pair[0] for pair in pairs]
    values_b = [pair[1] for pair in pairs]
    absolute = [abs(a - b) for a, b in pairs]
    relative = [abs(a - b) / abs(a) for a, b in pairs if a != 0]
    changed = sum(
        not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12) for a, b in pairs
    )
    result = {
        "pairs": len(pairs),
        "changed": changed,
        "pearson": safe_round(pearson(values_a, values_b)) if pairs else None,
        "spearman": (
            safe_round(pearson(average_ranks(values_a), average_ranks(values_b)))
            if pairs
            else None
        ),
        "median_absolute_delta": (
            safe_round(statistics.median(absolute)) if absolute else None
        ),
        "max_absolute_delta": safe_round(max(absolute)) if absolute else None,
        "median_relative_delta": (
            safe_round(statistics.median(relative)) if relative else None
        ),
        "missing_pairs": missing_both,
        "nonfinite_pairs": nonfinite_pairs,
        "nonfinite_value_mismatches": nonfinite_value_mismatches,
        "text_pairs_in_numeric_column": text_pairs,
        "text_value_mismatches": text_value_mismatches,
        "missing_mismatches": missing_mismatches,
        "nonfinite_mismatches": nonfinite_mismatches,
        "type_mismatches": type_mismatches,
        "mismatch_examples": mismatch_examples,
    }
    normalized_name = normalize_header(new_column)
    if pairs and ("log2fc" in normalized_name or "fold_change" in normalized_name):
        result["sign_agreement"] = safe_round(
            sum((a > 0) == (b > 0) and (a < 0) == (b < 0) for a, b in pairs)
            / len(pairs)
        )
    return result


def classify_numeric_cell(value: object) -> tuple[str, object]:
    if missing_value(value):
        return "missing", None
    number = numeric_value(value)
    if number is not None:
        return "finite", number
    if nonfinite_numeric_token(value):
        text = display_value(value).lower()
        if "nan" in text:
            return "nonfinite", "nan"
        return "nonfinite", "-infinity" if text.startswith("-") else "infinity"
    return "text", display_value(value)


def numeric_comparison_changed(comparison: dict) -> bool:
    return any(
        comparison[name]
        for name in (
            "changed",
            "nonfinite_value_mismatches",
            "text_value_mismatches",
            "missing_mismatches",
            "nonfinite_mismatches",
            "type_mismatches",
        )
    )


def compare_tables(old: dict, new: dict) -> dict:
    old_rows = old["keyed_rows"]
    new_rows = new["keyed_rows"]
    old_keys = set(old_rows)
    new_keys = set(new_rows)
    common_keys = sorted(old_keys & new_keys)
    union_size = len(old_keys | new_keys)

    old_columns = set(old["columns"])
    new_columns = set(new["columns"])
    common_columns = sorted(old_columns & new_columns)
    numeric = {}
    numeric_candidates = set(old["summary"]["numeric_columns"]) | set(
        new["summary"]["numeric_columns"]
    )
    for column in common_columns:
        if column in numeric_candidates and column not in set(old["key_columns"]) | set(
            new["key_columns"]
        ):
            numeric[column] = compare_numeric_column(
                old_rows, new_rows, common_keys, column, column
            )

    text_changed = 0
    compared_text = 0
    numeric_columns = set(numeric)
    for key in common_keys:
        for column in common_columns:
            if column in numeric_columns or column in set(old["key_columns"]) | set(
                new["key_columns"]
            ):
                continue
            compared_text += 1
            if display_value(old_rows[key].get(column)) != display_value(
                new_rows[key].get(column)
            ):
                text_changed += 1

    result = {
        "baseline_rows": len(old_rows),
        "candidate_rows": len(new_rows),
        "baseline_key_columns": list(old["key_columns"]),
        "candidate_key_columns": list(new["key_columns"]),
        "key_columns_match": old["key_columns"] == new["key_columns"],
        "common_keys": len(common_keys),
        "key_jaccard": safe_round(len(common_keys) / union_size) if union_size else 1.0,
        "baseline_only_keys": len(old_keys - new_keys),
        "baseline_only_key_examples": sorted(old_keys - new_keys)[:10],
        "candidate_only_keys": len(new_keys - old_keys),
        "candidate_only_key_examples": sorted(new_keys - old_keys)[:10],
        "baseline_only_columns": sorted(old_columns - new_columns),
        "candidate_only_columns": sorted(new_columns - old_columns),
        "text_cells_compared": compared_text,
        "text_cells_changed": text_changed,
        "numeric_columns": numeric,
    }
    result["semantic_equal"] = bool(
        result["key_columns_match"]
        and not result["baseline_only_keys"]
        and not result["candidate_only_keys"]
        and not result["baseline_only_columns"]
        and not result["candidate_only_columns"]
        and not result["text_cells_changed"]
        and not any(numeric_comparison_changed(item) for item in numeric.values())
    )
    return result


def public_table_summary(table: dict) -> dict:
    return table["summary"]


def compare_tabular_artifacts(old: dict, new: dict) -> dict:
    old_tables = old["tables"]
    new_tables = new["tables"]
    comparisons = {}
    for name in sorted(set(old_tables) | set(new_tables)):
        if name not in old_tables:
            comparisons[name] = {
                "status": "candidate_only",
                "candidate": public_table_summary(new_tables[name]),
            }
        elif name not in new_tables:
            comparisons[name] = {
                "status": "baseline_only",
                "baseline": public_table_summary(old_tables[name]),
            }
        else:
            detail = compare_tables(old_tables[name], new_tables[name])
            comparisons[name] = {
                "status": "unchanged" if detail["semantic_equal"] else "changed",
                "baseline": public_table_summary(old_tables[name]),
                "candidate": public_table_summary(new_tables[name]),
                "comparison": detail,
            }
    return {
        "semantic_equal": set(old_tables) == set(new_tables)
        and all(
            result["status"] == "unchanged" for result in comparisons.values()
        ),
        "tables": comparisons,
    }


def compare_gff(old: dict, new: dict) -> dict:
    old_counter = Counter(old["keys"])
    new_counter = Counter(new["keys"])
    old_keys = set(old_counter)
    new_keys = set(new_counter)
    common = sum((old_counter & new_counter).values())
    union = sum((old_counter | new_counter).values())
    tolerant = tolerant_gff_matches(old["records"], new["records"])
    return {
        "semantic_equal": Counter(map(gff_semantic_key, old["records"]))
        == Counter(map(gff_semantic_key, new["records"])),
        "baseline": old["summary"],
        "candidate": new["summary"],
        "exact_coordinate_matches": common,
        "exact_coordinate_jaccard": safe_round(common / union) if union else 1.0,
        "within_3nt_matches": tolerant,
        "within_3nt_baseline_fraction": (
            safe_round(tolerant / len(old["records"])) if old["records"] else 1.0
        ),
        "within_3nt_candidate_fraction": (
            safe_round(tolerant / len(new["records"])) if new["records"] else 1.0
        ),
        "baseline_only_coordinate_examples": [
            list(item) for item in sorted(old_keys - new_keys)[:10]
        ],
        "candidate_only_coordinate_examples": [
            list(item) for item in sorted(new_keys - old_keys)[:10]
        ],
    }


def json_semantically_equal(old: object, new: object) -> bool:
    """Compare decoded JSON without conflating booleans and numbers."""

    if isinstance(old, bool) or isinstance(new, bool):
        return isinstance(old, bool) and isinstance(new, bool) and old is new
    if isinstance(old, (int, float)) and isinstance(new, (int, float)):
        # JSON has one number type, so harmless integer/decimal spelling
        # differences remain equivalent. Non-finite values are rejected while
        # parsing and booleans were handled separately above.
        return old == new
    if type(old) is not type(new):
        return False
    if isinstance(old, dict):
        return old.keys() == new.keys() and all(
            json_semantically_equal(old[key], new[key]) for key in old
        )
    if isinstance(old, list):
        return len(old) == len(new) and all(
            json_semantically_equal(old_item, new_item)
            for old_item, new_item in zip(old, new, strict=True)
        )
    return old == new


def compare_json(old: dict, new: dict) -> dict:
    old_value = old["value"]
    new_value = new["value"]
    equal = json_semantically_equal(old_value, new_value)
    return {
        "semantic_equal": equal,
        "baseline": old["summary"],
        "candidate": new["summary"],
        "equal": equal,
        "baseline_only_top_level_keys": (
            sorted(set(old_value) - set(new_value))
            if isinstance(old_value, dict) and isinstance(new_value, dict)
            else []
        ),
        "candidate_only_top_level_keys": (
            sorted(set(new_value) - set(old_value))
            if isinstance(old_value, dict) and isinstance(new_value, dict)
            else []
        ),
    }


def compare_bam(old: dict, new: dict) -> dict:
    comparison = compare_tabular_artifacts(old, new)
    signature_equal = old["semantic_signature"] == new["semantic_signature"]
    comparison["header_equal"] = (
        old["semantic_signature"]["stable_header_sha256"]
        == new["semantic_signature"]["stable_header_sha256"]
    )
    comparison["program_provenance_equal"] = (
        old["provenance_sha256"] == new["provenance_sha256"]
    )
    comparison["alignment_digest_equal"] = (
        old["semantic_signature"]["alignment_digest"]
        == new["semantic_signature"]["alignment_digest"]
    )
    comparison["semantic_equal"] = comparison["semantic_equal"] and signature_equal
    return comparison


def compare_bigwig(old: dict, new: dict) -> dict:
    def decoded_semantics(summary: dict) -> dict:
        return {
            "field_count": summary["field_count"],
            "defined_field_count": summary["defined_field_count"],
            "zoom_reductions": [item["reduction"] for item in summary["zoom_levels"]],
            "chromosomes": summary["chromosomes"],
            "total_summary": summary["total_summary"],
        }

    summary_equal = decoded_semantics(old["summary"]) == decoded_semantics(
        new["summary"]
    )
    byte_digest_equal = old["sha256"] == new["sha256"]
    return {
        "semantic_equal": byte_digest_equal,
        "comparison_mode": (
            "conservative compressed-byte equality with decoded BigWig metadata"
        ),
        "decoded_summary_equal": summary_equal,
        "byte_digest_equal": byte_digest_equal,
        "baseline": old["summary"],
        "candidate": new["summary"],
    }


def comparison_hint(path: str) -> str | None:
    lower = path.lower()
    if "metageneprofiling" in lower and "global_readcounts" in lower:
        return (
            "Global metagene slices intentionally changed: the old release dropped "
            "minus/start and plus/stop windows."
        )
    if "metageneprofiling" in lower and "readcounts_stop" in lower:
        return (
            "Stop-codon profiles intentionally changed orientation so both strands "
            "run from coding sequence toward downstream sequence."
        )
    if lower.startswith(("xtail/", "riborex/", "deltate/")):
        return (
            "Differential results may change because engines, contrast handling and "
            "failure validation were corrected; review direction and correlation."
        )
    if "prediction" in lower or "deepribo" in lower or "reparation" in lower:
        return (
            "Prediction results may change after parser, strand geometry, aggregation "
            "and engine-boundary fixes; review exact and within-3-nt overlap."
        )
    return None


def is_allowed_missing(path: str, patterns: tuple[str, ...]) -> bool:
    return any(fnmatch.fnmatch(path, pattern) for pattern in patterns)


def add_issue(report: dict, severity: str, path: str, message: str) -> None:
    report["issues"].append(
        {"severity": severity, "path": path, "message": message}
    )


def validate_thresholds(
    report: dict,
    path: str,
    comparison: dict,
    minimum_key_overlap: float | None,
    minimum_correlation: float | None,
) -> None:
    if minimum_key_overlap is None and minimum_correlation is None:
        return
    for sheet, sheet_result in comparison.get("tables", {}).items():
        detail = sheet_result.get("comparison", {})
        overlap = detail.get("key_jaccard")
        if (
            minimum_key_overlap is not None
            and overlap is not None
            and overlap < minimum_key_overlap
        ):
            add_issue(
                report,
                "error",
                path,
                f"sheet {sheet!r} key Jaccard {overlap} is below "
                f"{minimum_key_overlap}",
            )
        if minimum_correlation is None:
            continue
        for column, metric in detail.get("numeric_columns", {}).items():
            correlation = metric.get("spearman")
            if metric.get("pairs", 0) >= 3 and correlation is None:
                # A constant column has no defined rank correlation even when
                # every finite value is identical between the two runs.  Such
                # a column carries no drift to threshold; changed degenerate
                # columns must still fail because their undefined correlation
                # cannot demonstrate the requested agreement.
                if not numeric_comparison_changed(metric):
                    continue
                add_issue(
                    report,
                    "error",
                    path,
                    f"sheet {sheet!r} column {column!r} Spearman is undefined "
                    f"for {metric['pairs']} finite pairs",
                )
                continue
            if (
                correlation is not None
                and metric.get("pairs", 0) >= 3
                and correlation < minimum_correlation
            ):
                add_issue(
                    report,
                    "error",
                    path,
                    f"sheet {sheet!r} column {column!r} Spearman {correlation} "
                    f"is below {minimum_correlation}",
                )


def compare_runs(
    baseline_root: Path,
    candidate_root: Path,
    *,
    allow_missing: tuple[str, ...] = (),
    minimum_key_overlap: float | None = None,
    minimum_correlation: float | None = None,
) -> dict:
    """Inspect and compare two result roots, returning a JSON-safe report."""

    for name, value in (
        ("minimum_key_overlap", minimum_key_overlap),
        ("minimum_correlation", minimum_correlation),
    ):
        if value is not None and not 0 <= value <= 1:
            raise ComparisonError(f"{name} must be between 0 and 1")

    baseline_root = resolve_path(
        baseline_root, description=f"baseline run root {baseline_root}"
    )
    candidate_root = resolve_path(
        candidate_root, description=f"candidate run root {candidate_root}"
    )
    if baseline_root == candidate_root:
        raise ComparisonError("baseline and candidate roots resolve to the same directory")
    baseline_paths = discover_outputs(baseline_root)
    candidate_paths = discover_outputs(candidate_root)
    if not baseline_paths:
        raise ComparisonError(f"no supported outputs found under {baseline_root}")
    if not candidate_paths:
        raise ComparisonError(f"no supported outputs found under {candidate_root}")

    report = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "baseline_root": str(baseline_root),
        "candidate_root": str(candidate_root),
        "settings": {
            "allow_missing": list(allow_missing),
            "minimum_key_overlap": minimum_key_overlap,
            "minimum_correlation": minimum_correlation,
        },
        "summary": {},
        "issues": [],
        "artifacts": [],
    }

    exact = changed = candidate_only = baseline_only = invalid = 0
    all_paths = sorted(set(baseline_paths) | set(candidate_paths))
    for relative in all_paths:
        entry = {"path": relative}
        old_path = baseline_paths.get(relative)
        new_path = candidate_paths.get(relative)
        if old_path is None:
            candidate_only += 1
            entry["status"] = "candidate_only"
            try:
                inspected = inspect_artifact(new_path, relative)
                entry["kind"] = inspected["kind"]
                entry["candidate"] = artifact_public_summary(inspected)
            except (ComparisonError, OSError) as error:
                invalid += 1
                entry["status"] = "invalid"
                entry["error"] = str(error)
                add_issue(report, "error", relative, f"invalid candidate: {error}")
            report["artifacts"].append(entry)
            continue
        if new_path is None:
            baseline_only += 1
            entry["status"] = "baseline_only"
            entry["allowed"] = is_allowed_missing(relative, allow_missing)
            if not entry["allowed"]:
                add_issue(report, "error", relative, "candidate is missing baseline output")
            report["artifacts"].append(entry)
            continue

        try:
            old = inspect_artifact(
                old_path, relative, baseline_compatibility=True
            )
            new = inspect_artifact(new_path, relative)
        except (ComparisonError, OSError) as error:
            invalid += 1
            entry["status"] = "invalid"
            entry["error"] = str(error)
            add_issue(report, "error", relative, str(error))
            report["artifacts"].append(entry)
            continue
        if old["kind"] != new["kind"]:
            invalid += 1
            entry["status"] = "invalid"
            entry["error"] = f"type changed from {old['kind']} to {new['kind']}"
            add_issue(report, "error", relative, entry["error"])
            report["artifacts"].append(entry)
            continue

        byte_identical = (
            old["sha256"] == new["sha256"]
            if "sha256" in old and "sha256" in new
            else None
        )
        entry["kind"] = old["kind"]
        entry["baseline"] = artifact_public_summary(old)
        entry["candidate"] = artifact_public_summary(new)
        entry["byte_identical"] = byte_identical

        legacy_compatibility = old.get("summary", {}).get(
            "legacy_deepribo_phase_compatibility"
        )
        if legacy_compatibility is not None:
            entry["baseline_compatibility"] = legacy_compatibility
            add_issue(
                report,
                "warning",
                relative,
                "normalized known HRIBO 1.8 DeepRibo CDS phase misuse to phase "
                f"0 for {legacy_compatibility['records_normalized']} baseline "
                "record(s); original values remain in artifact metadata",
            )
        legacy_empty_compatibility = old.get("summary", {}).get(
            "legacy_zero_byte_gff_compatibility"
        )
        if legacy_empty_compatibility is not None:
            entry["baseline_zero_result_compatibility"] = (
                legacy_empty_compatibility
            )
            add_issue(
                report,
                "warning",
                relative,
                "interpreted a known HRIBO 1.8 zero-byte per-condition GFF "
                "as an empty baseline feature set",
            )

        if old["kind"] in {"xlsx", "csv", "tsv", "tab"}:
            entry["comparison"] = compare_tabular_artifacts(old, new)
        elif old["kind"] == "bam":
            entry["comparison"] = compare_bam(old, new)
        elif old["kind"] == "bigwig":
            entry["comparison"] = compare_bigwig(old, new)
        elif old["kind"] == "gff":
            entry["comparison"] = compare_gff(old, new)
        elif old["kind"] == "json":
            entry["comparison"] = compare_json(old, new)

        semantic_equal = entry["comparison"]["semantic_equal"]
        entry["status"] = "unchanged" if semantic_equal else "changed"
        if semantic_equal:
            exact += 1
        else:
            changed += 1
            hint = comparison_hint(relative)
            if hint:
                entry["review_hint"] = hint

        validate_thresholds(
            report,
            relative,
            entry.get("comparison", {}),
            minimum_key_overlap,
            minimum_correlation,
        )
        report["artifacts"].append(entry)

    compared = exact + changed
    if compared == 0 and invalid == 0:
        add_issue(report, "error", "", "the runs have no comparable output paths")
    error_count = sum(issue["severity"] == "error" for issue in report["issues"])
    warning_count = sum(
        issue["severity"] == "warning" for issue in report["issues"]
    )
    report["summary"] = {
        "baseline_artifacts": len(baseline_paths),
        "candidate_artifacts": len(candidate_paths),
        "comparable_artifacts": compared,
        "unchanged_artifacts": exact,
        "changed_artifacts": changed,
        "candidate_only_artifacts": candidate_only,
        "baseline_only_artifacts": baseline_only,
        "invalid_artifacts": invalid,
        "errors": error_count,
        "warnings": warning_count,
        "review_required": bool(
            changed or candidate_only or baseline_only or invalid or warning_count
        ),
    }
    return report


def artifact_public_summary(artifact: dict) -> dict:
    result = {"bytes": artifact["bytes"]}
    if "sha256" in artifact:
        result["sha256"] = artifact["sha256"]
    if "summary" in artifact:
        result.update(artifact["summary"])
    if "tables" in artifact:
        result["tables"] = {
            name: public_table_summary(table)
            for name, table in artifact["tables"].items()
        }
    return result


def write_json_atomic(path: Path, report: dict) -> None:
    path = resolve_path(path, description=f"report path {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def parse_arguments(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare primary scientific outputs from a baseline and candidate "
            "HRIBO result directory. Changed values are reported for review; "
            "missing or malformed candidate artifacts fail by default."
        )
    )
    parser.add_argument("baseline", type=Path, help="HRIBO baseline result directory")
    parser.add_argument("candidate", type=Path, help="HRIBO candidate result directory")
    parser.add_argument(
        "--report", type=Path, required=True, help="destination JSON report"
    )
    parser.add_argument(
        "--allow-missing",
        action="append",
        default=[],
        metavar="GLOB",
        help="allow an intentional candidate-missing relative path (repeatable)",
    )
    parser.add_argument(
        "--minimum-key-overlap",
        type=float,
        help="fail when a comparable table's key Jaccard is below this value",
    )
    parser.add_argument(
        "--minimum-correlation",
        type=float,
        help=(
            "fail when a numeric column with at least three pairs has a Spearman "
            "correlation below this value"
        ),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_arguments(argv)
    try:
        report = compare_runs(
            args.baseline,
            args.candidate,
            allow_missing=tuple(args.allow_missing),
            minimum_key_overlap=args.minimum_key_overlap,
            minimum_correlation=args.minimum_correlation,
        )
        write_json_atomic(args.report, report)
    except (ComparisonError, OSError) as error:
        print(f"Comparison failed: {error}", file=sys.stderr)
        return 2

    summary = report["summary"]
    print(
        "Compared "
        f"{summary['comparable_artifacts']} artifacts: "
        f"{summary['unchanged_artifacts']} unchanged, "
        f"{summary['changed_artifacts']} changed, "
        f"{summary['candidate_only_artifacts']} candidate-only, "
        f"{summary['baseline_only_artifacts']} baseline-only, "
        f"{summary['invalid_artifacts']} invalid."
    )
    print(f"Report: {args.report.resolve()}")
    if summary["errors"]:
        print(f"Comparison has {summary['errors']} structural/threshold error(s).")
        return 1
    if summary["review_required"]:
        print("Scientific differences are present and require review.")
    else:
        print("All discovered outputs are unchanged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
