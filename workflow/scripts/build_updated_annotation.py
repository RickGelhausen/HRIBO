#!/usr/bin/env python3
"""Build one valid GFF3 annotation from reference and prediction records.

Reference attributes are copied byte-for-byte. Prediction IDs are prefixed
with their record source and allocated around reference/encoded-name
collisions; prediction-internal Parent and Derives_from references are updated
to match.
"""

import argparse
import string
from dataclasses import dataclass
from pathlib import Path

from concatenate_gff import atomic_write, render, sort_records


REFERENCE_ATTRIBUTES = {"Parent", "Derives_from"}
GFF_ID_SAFE = frozenset(string.ascii_letters + string.digits + ".^*$@!+_?-|")
HEX_DIGITS = frozenset(string.hexdigits)


@dataclass
class Record:
    fields: list[str]
    start: int
    end: int
    path: Path
    line_number: int

    @property
    def source(self):
        return self.fields[1]

    @property
    def context(self):
        return f"{self.path}:{self.line_number}"

    def sortable(self):
        return self.fields, self.start, self.end


def _coordinate(value, path, line_number, column):
    try:
        coordinate = int(value)
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number}: {column} coordinate must be an integer, "
            f"found {value!r}"
        ) from error
    if coordinate < 1:
        raise ValueError(
            f"{path}:{line_number}: {column} coordinate must be positive, "
            f"found {coordinate}"
        )
    return coordinate


def read_gff3(path):
    """Read feature records while rejecting malformed structural data."""

    path = Path(path)
    records = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                if stripped.startswith("##gff-version") and stripped != "##gff-version 3":
                    raise ValueError(
                        f"{path}:{line_number}: expected '##gff-version 3', "
                        f"found {stripped!r}"
                    )
                if stripped == "##FASTA":
                    raise ValueError(
                        f"{path}:{line_number}: embedded FASTA is not supported"
                    )
                continue

            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{path}:{line_number}: expected 9 tab-separated GFF3 columns, "
                    f"found {len(fields)}"
                )
            if not fields[0] or not fields[1] or not fields[2]:
                raise ValueError(
                    f"{path}:{line_number}: seqid, source, and type must be non-empty"
                )

            start = _coordinate(fields[3], path, line_number, "start")
            end = _coordinate(fields[4], path, line_number, "end")
            if start > end:
                raise ValueError(
                    f"{path}:{line_number}: start coordinate {start} exceeds end {end}"
                )
            if fields[6] not in {"+", "-", ".", "?"}:
                raise ValueError(
                    f"{path}:{line_number}: invalid strand {fields[6]!r}"
                )
            if fields[7] not in {"0", "1", "2", "."}:
                raise ValueError(
                    f"{path}:{line_number}: invalid phase {fields[7]!r}"
                )

            records.append(Record(fields, start, end, path, line_number))
    return records


def parse_attributes(record):
    """Return ordered attributes and whether the source used a final semicolon."""

    attributes = record.fields[8]
    if attributes == ".":
        return [], False
    if not attributes:
        raise ValueError(f"{record.context}: attribute column must not be empty")

    trailing_semicolon = attributes.endswith(";")
    fields = attributes.split(";")
    if trailing_semicolon:
        fields.pop()
    if any(not field for field in fields):
        raise ValueError(f"{record.context}: empty field in attribute column")

    pairs = []
    for field in fields:
        if "=" not in field:
            raise ValueError(
                f"{record.context}: malformed GFF3 attribute {field!r}; expected key=value"
            )
        key, value = field.split("=", 1)
        if not key:
            raise ValueError(f"{record.context}: GFF3 attribute key must not be empty")
        pairs.append((key, value))
    return pairs, trailing_semicolon


def _escape_id_component(value):
    """Percent-encode a namespace component without double-encoding GFF3 escapes."""

    encoded = []
    index = 0
    while index < len(value):
        character = value[index]
        if (
            character == "%"
            and index + 2 < len(value)
            and value[index + 1] in HEX_DIGITS
            and value[index + 2] in HEX_DIGITS
        ):
            encoded.append(value[index : index + 3])
            index += 3
            continue

        for byte in character.encode("utf-8"):
            byte_character = chr(byte)
            if byte_character in GFF_ID_SAFE:
                encoded.append(byte_character)
            else:
                encoded.append(f"%{byte:02X}")
        index += 1
    return "".join(encoded)


def namespaced_id(source, identifier):
    if not identifier:
        raise ValueError("prediction ID must not be empty")
    return f"{_escape_id_component(source)}:{_escape_id_component(identifier)}"


def _id_pair(pairs, record, record_kind):
    identifiers = [value for key, value in pairs if key == "ID"]
    if len(identifiers) > 1:
        raise ValueError(f"{record.context}: record contains more than one ID attribute")
    if identifiers and not identifiers[0]:
        raise ValueError(f"{record.context}: {record_kind} ID must not be empty")
    if identifiers and "," in identifiers[0]:
        raise ValueError(
            f"{record.context}: {record_kind} ID must contain one value"
        )
    return identifiers[0] if identifiers else None


def _unique_namespaced_ids(keys, reserved_ids):
    """Allocate deterministic IDs that cannot collide after percent encoding."""

    allocated = {}
    used = set(reserved_ids)
    for source, identifier in sorted(keys):
        base = namespaced_id(source, identifier)
        candidate = base
        suffix = 1
        while candidate in used:
            label = "prediction" if suffix == 1 else f"prediction{suffix}"
            candidate = f"{base}:{label}"
            suffix += 1
        allocated[(source, identifier)] = candidate
        used.add(candidate)
    return allocated


def _definitions(predictions, reserved_ids):
    """Index prediction IDs by identifier/source with globally unique values."""

    parsed = []
    keys = set()
    for record in predictions:
        pairs, trailing_semicolon = parse_attributes(record)
        identifier = _id_pair(pairs, record, "prediction")
        parsed.append((record, pairs, trailing_semicolon, identifier))
        if identifier is not None:
            keys.add((record.source, identifier))

    allocated = _unique_namespaced_ids(keys, reserved_ids)
    definitions = {}
    for (source, identifier), updated_identifier in allocated.items():
        definitions.setdefault(identifier, {})[source] = updated_identifier
    return definitions, parsed


def _rewrite_reference(value, record, definitions, attribute):
    if not value:
        raise ValueError(
            f"{record.context}: {attribute} reference must not be empty"
        )

    rewritten = []
    for identifier in value.split(","):
        if not identifier:
            raise ValueError(
                f"{record.context}: {attribute} contains an empty reference"
            )
        candidates = definitions.get(identifier)
        if not candidates:
            rewritten.append(identifier)
        elif record.source in candidates:
            rewritten.append(candidates[record.source])
        elif len(candidates) == 1:
            rewritten.append(next(iter(candidates.values())))
        else:
            sources = ", ".join(sorted(candidates))
            raise ValueError(
                f"{record.context}: ambiguous internal {attribute} reference "
                f"{identifier!r}; it is defined by sources {sources}"
            )
    return ",".join(rewritten)


def namespace_predictions(predictions, reserved_ids=()):
    definitions, parsed = _definitions(predictions, reserved_ids)
    updated = []
    for record, pairs, trailing_semicolon, identifier in parsed:
        rewritten_pairs = []
        for key, value in pairs:
            if key == "ID":
                value = definitions[identifier][record.source]
            elif key in REFERENCE_ATTRIBUTES:
                value = _rewrite_reference(value, record, definitions, key)
            rewritten_pairs.append((key, value))

        if rewritten_pairs:
            attributes = ";".join(
                f"{key}={value}" for key, value in rewritten_pairs
            )
            if trailing_semicolon:
                attributes += ";"
        else:
            attributes = "."
        fields = [*record.fields[:8], attributes]
        updated.append(Record(fields, record.start, record.end, record.path, record.line_number))
    return updated


def record_ids(records, record_kind):
    """Collect IDs without changing the original attribute serialization."""

    identifiers = set()
    for record in records:
        pairs, _trailing_semicolon = parse_attributes(record)
        identifier = _id_pair(pairs, record, record_kind)
        if identifier is not None:
            identifiers.add(identifier)
    return identifiers


def build_updated_annotation(annotation, predictions, output):
    reference_records = read_gff3(annotation)
    prediction_records = []
    for prediction in predictions:
        prediction_records.extend(read_gff3(prediction))

    updated_predictions = namespace_predictions(
        prediction_records, record_ids(reference_records, "reference")
    )
    records = [record.sortable() for record in reference_records + updated_predictions]
    atomic_write(output, render(sort_records(records)))


def main():
    parser = argparse.ArgumentParser(
        description="Combine a checked annotation with namespaced prediction GFF3 records."
    )
    parser.add_argument("--annotation", required=True, help="checked reference GFF3")
    parser.add_argument(
        "--predictions",
        nargs="*",
        default=[],
        metavar="GFF3",
        help="zero or more prediction GFF3 files",
    )
    parser.add_argument("--output", required=True, help="combined output GFF3")
    args = parser.parse_args()

    try:
        build_updated_annotation(args.annotation, args.predictions, args.output)
    except (OSError, ValueError) as error:
        parser.exit(1, f"build_updated_annotation.py: error: {error}\n")


if __name__ == "__main__":
    main()
