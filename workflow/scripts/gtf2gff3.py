#!/usr/bin/env python
"""Convert a consistently formatted GTF annotation to GFF3.

GFF3 input is copied unchanged. GTF records are grouped by ``gene_id`` and
receive deterministic GFF3 IDs while retaining their remaining attributes.
The converter deliberately rejects mixed or unidentifiable attribute syntax:
silently guessing here would corrupt every annotation-derived result later in
the workflow.
"""

from __future__ import annotations

import argparse
import os
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from shutil import copyfile

import gff_utils


RNA_FEATURES = {"ncrna", "rrna", "trna", "srna"}
IGNORED_GTF_FEATURES = {"exon", "start_codon", "stop_codon", "transcript"}
GENERATED_ATTRIBUTE_KEYS = {"gene_id", "id", "locus_tag", "parent"}
GFF3_RESERVED_ATTRIBUTE_CHARACTERS = frozenset("%;=&,")


class AnnotationConversionError(ValueError):
    """An annotation cannot be converted without guessing its meaning."""


@dataclass(frozen=True)
class Record:
    seq_name: str
    source: str
    feature: str
    start: int
    stop: int
    score: str
    strand: str
    phase: str
    attributes: str
    line_number: int

    def to_gff_line(self) -> str:
        return "\t".join(
            (
                self.seq_name,
                self.source,
                self.feature,
                str(self.start),
                str(self.stop),
                self.score,
                self.strand,
                self.phase,
                self.attributes,
            )
        )


def read_annotation(path: Path) -> tuple[list[Record], bool]:
    """Read the nine GFF/GTF columns without coercing identifiers or dots."""
    records = []
    declares_gff3 = False
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\n\r")
            if not stripped:
                continue
            if stripped.startswith("#"):
                declares_gff3 = declares_gff3 or stripped == "##gff-version 3"
                continue

            fields = stripped.split("\t")
            if len(fields) != 9:
                raise AnnotationConversionError(
                    f"line {line_number} has {len(fields)} columns; expected nine tab-separated columns"
                )
            try:
                start, stop = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise AnnotationConversionError(
                    f"line {line_number} has non-integer coordinates {fields[3]!r}-{fields[4]!r}"
                ) from exc
            if start < 1 or stop < start:
                raise AnnotationConversionError(
                    f"line {line_number} has invalid coordinates {start}-{stop}"
                )

            records.append(
                Record(
                    fields[0],
                    fields[1],
                    fields[2],
                    start,
                    stop,
                    fields[5],
                    fields[6],
                    fields[7],
                    fields[8],
                    line_number,
                )
            )

    if not records:
        raise AnnotationConversionError(f"{path} contains no annotation records")
    return records, declares_gff3


def split_attribute_fields(attributes: str) -> list[str]:
    """Split column nine on semicolons that are outside GTF quoted values."""
    fields = []
    current = []
    quoted = False
    for character in attributes:
        if character == '"':
            quoted = not quoted
        if character == ";" and not quoted:
            field = "".join(current).strip()
            if field:
                fields.append(field)
            current = []
        else:
            current.append(character)

    field = "".join(current).strip()
    if field:
        fields.append(field)
    return fields


def split_gtf_attributes(attributes: str) -> list[tuple[str, str]]:
    """Return already-validated GTF attributes without splitting quoted text."""
    pairs = []
    for field in split_attribute_fields(attributes):
        match = gff_utils.GTF2_PAIR.match(field)
        if match is not None:
            pairs.append((match.group("key"), match.group("value")))
    return pairs


def parse_gtf_attributes(attributes: str) -> dict[str, str]:
    """Parse GTF attributes case-insensitively, retaining the first value."""
    parsed = {}
    for key, value in split_gtf_attributes(attributes):
        parsed.setdefault(key.lower(), value)
    return parsed


def attribute_syntax(attributes: str, line_number: int) -> set[str]:
    """Return the syntaxes used by one attribute column."""
    if attributes.strip() in {"", "."}:
        return set()

    syntaxes = set()
    for field in split_attribute_fields(attributes):
        if gff_utils.GTF2_PAIR.match(field):
            syntaxes.add("GTF")
        elif "=" in field and field.split("=", 1)[0].strip():
            syntaxes.add("GFF3")
        else:
            raise AnnotationConversionError(
                f"line {line_number} has unsupported attribute {field!r}; "
                'use GTF key "value" or GFF3 key=value syntax'
            )
    return syntaxes


def annotation_format(records: list[Record], declares_gff3: bool) -> str:
    """Identify one consistent format and require its defining identifier."""
    syntaxes = set()
    for record in records:
        syntaxes.update(attribute_syntax(record.attributes, record.line_number))

    if len(syntaxes) > 1:
        raise AnnotationConversionError(
            "mixed GFF3 and GTF attribute syntax; use one format consistently"
        )
    if not syntaxes:
        raise AnnotationConversionError(
            'no supported identifiers found; expected GTF gene_id "..." or GFF3 ID=...'
        )

    detected = next(iter(syntaxes))
    if declares_gff3 and detected != "GFF3":
        raise AnnotationConversionError(
            "the file declares GFF3 but its records use GTF attribute syntax"
        )

    parser = parse_gtf_attributes if detected == "GTF" else gff_utils.parse_attributes
    parsed = [parser(record.attributes) for record in records]
    if detected == "GFF3" and not any("id" in attributes for attributes in parsed):
        raise AnnotationConversionError(
            "GFF3 input contains no ID attribute; HRIBO cannot identify its features"
        )
    if detected == "GTF" and not any("gene_id" in attributes for attributes in parsed):
        raise AnnotationConversionError(
            'GTF input contains no gene_id attribute; expected gene_id "..."'
        )
    return detected


def converted_attributes(
    attributes: str,
    identifier: str,
    gene_id: str,
    locus_tag: str,
    parent: str = "",
) -> str:
    """Add one canonical identity set while retaining lowercase GTF provenance."""
    remaining = [
        (escape_gff3_component(key), escape_gff3_component(value))
        for key, value in split_gtf_attributes(attributes)
        if key.lower() not in GENERATED_ATTRIBUTE_KEYS
    ]
    prefix = [
        ("ID", escape_gff3_component(identifier)),
        ("locus_tag", escape_gff3_component(locus_tag)),
    ]
    if parent:
        prefix.append(("Parent", escape_gff3_component(parent)))
    prefix.append(("gene_id", escape_gff3_component(gene_id)))
    return gff_utils.format_attributes(prefix + remaining)


def escape_gff3_component(value: str) -> str:
    """Percent-encode one GTF attribute component for GFF3 column nine.

    Column nine reserves ``;``, ``=``, ``&`` and ``,`` as separators and uses
    URL-style escapes for those characters, literal percent signs, ASCII
    controls and DEL.  Spaces, punctuation and UTF-8 are valid unescaped text.
    A percent triplet in GTF is literal source text, so its percent sign is
    encoded too; already-converted GFF3 input takes the separate pass-through
    path and never reaches this function.
    """
    escaped = []
    for character in value:
        if (
            character in GFF3_RESERVED_ATTRIBUTE_CHARACTERS
            or ord(character) < 32
            or ord(character) == 127
        ):
            escaped.extend(f"%{byte:02X}" for byte in character.encode("utf-8"))
        else:
            escaped.append(character)
    return "".join(escaped)


def structural_record_key(record: Record) -> tuple:
    """Stable biological/file-content order, independent of input row order."""
    return (
        record.seq_name,
        record.start,
        record.stop,
        record.strand,
        record.feature.lower(),
        record.feature,
        record.source,
        record.score,
        record.phase,
        record.attributes,
    )


def output_record_key(record: Record) -> tuple:
    """Keep parents before same-start children and resolve every remaining tie."""
    parent_rank = 0 if record.feature.lower() in {"gene", "pseudogene"} else 1
    return (
        record.seq_name,
        record.start,
        parent_rank,
        record.stop,
        record.strand,
        record.feature.lower(),
        record.feature,
        record.source,
        record.score,
        record.phase,
        record.attributes,
    )


def group_locus_tags(records_by_gene: dict[str, list[Record]]) -> dict[str, str]:
    """Choose one explicit nonempty locus tag per gene, or use ``gene_id``."""
    resolved = {}
    for gene_id, records in records_by_gene.items():
        occurrences = defaultdict(list)
        for record in records:
            for key, value in split_gtf_attributes(record.attributes):
                if key.lower() == "locus_tag" and value:
                    occurrences[value].append(record.line_number)

        if len(occurrences) > 1:
            conflicts = ", ".join(
                f"{value!r} (lines {', '.join(map(str, occurrences[value]))})"
                for value in sorted(occurrences)
            )
            raise AnnotationConversionError(
                f"gene_id {gene_id!r} has conflicting locus_tag values: {conflicts}"
            )
        resolved[gene_id] = next(iter(occurrences), gene_id)
    return resolved


def validate_gene_groups(
    records_by_gene: dict[str, list[Record]], genes: dict[str, Record]
) -> None:
    """Reject groupings that cannot form a valid GFF3 parent hierarchy."""
    for gene_id, records in records_by_gene.items():
        seq_names = {record.seq_name for record in records}
        if len(seq_names) > 1:
            lines = ", ".join(str(record.line_number) for record in records)
            raise AnnotationConversionError(
                f"gene_id {gene_id!r} spans multiple sequence IDs on lines {lines}"
            )

        strands = {record.strand for record in records}
        if len(strands) > 1:
            lines = ", ".join(str(record.line_number) for record in records)
            raise AnnotationConversionError(
                f"gene_id {gene_id!r} spans multiple strands on lines {lines}"
            )

        gene = genes.get(gene_id)
        if gene is None:
            continue
        for record in records:
            if record is gene:
                continue
            if record.start < gene.start or record.stop > gene.stop:
                raise AnnotationConversionError(
                    f"line {record.line_number} ({record.feature} {record.start}-{record.stop}) "
                    f"lies outside gene_id {gene_id!r} on line {gene.line_number} "
                    f"({gene.start}-{gene.stop})"
                )


def group_gtf_records(records: list[Record]):
    genes = {}
    coding = defaultdict(list)
    rna = defaultdict(list)
    unknown = []
    records_by_gene = defaultdict(list)

    for record in records:
        attributes = parse_gtf_attributes(record.attributes)
        gene_id = attributes.get("gene_id", "")
        if not gene_id:
            raise AnnotationConversionError(
                f"line {record.line_number} has no gene_id; every convertible GTF record needs one"
            )
        records_by_gene[gene_id].append(record)

        feature = record.feature.lower()
        if feature in {"gene", "pseudogene"}:
            if gene_id in genes:
                raise AnnotationConversionError(
                    f"duplicate gene records for gene_id {gene_id!r} on lines "
                    f"{genes[gene_id].line_number} and {record.line_number}"
                )
            genes[gene_id] = record
        elif feature == "cds":
            coding[gene_id].append(record)
        elif feature in RNA_FEATURES:
            rna[gene_id].append(record)
        elif feature == "transcript":
            biotype = attributes.get("gene_biotype", "")
            if biotype.lower() in RNA_FEATURES:
                rna[gene_id].append(replace(record, feature=biotype))
        elif feature not in IGNORED_GTF_FEATURES:
            unknown.append((gene_id, record))

    validate_gene_groups(records_by_gene, genes)
    locus_tags = group_locus_tags(records_by_gene)
    return genes, coding, rna, unknown, locus_tags


def convert_gtf(records: list[Record]) -> list[Record]:
    genes, coding, rna, unknown, locus_tags = group_gtf_records(records)
    for children in coding.values():
        children.sort(key=structural_record_key)
    for children in rna.values():
        children.sort(key=structural_record_key)
    unknown.sort(key=lambda item: (structural_record_key(item[1]), item[0]))
    gene_ids = sorted(set(genes) | set(coding) | set(rna))

    converted = []
    gene_count = cds_count = rna_count = 1
    for gene_id in gene_ids:
        parent_id = f"gene{gene_count}"
        locus_tag = locus_tags[gene_id]
        if gene_id in genes:
            gene = genes[gene_id]
            converted.append(
                replace(
                    gene,
                    phase=".",
                    attributes=converted_attributes(
                        gene.attributes, parent_id, gene_id, locus_tag
                    ),
                )
            )

            for child in coding.get(gene_id, []):
                converted.append(
                    replace(
                        child,
                        feature="CDS",
                        attributes=converted_attributes(
                            child.attributes,
                            f"cds{cds_count}",
                            gene_id,
                            locus_tag,
                            parent_id,
                        ),
                    )
                )
                cds_count += 1

            for child in rna.get(gene_id, []):
                converted.append(
                    replace(
                        child,
                        phase=".",
                        attributes=converted_attributes(
                            child.attributes,
                            f"rna{rna_count}",
                            gene_id,
                            locus_tag,
                            parent_id,
                        ),
                    )
                )
                rna_count += 1

        else:
            coding_children = coding.get(gene_id, [])
            rna_children = rna.get(gene_id, [])
            children = coding_children + rna_children
            first = min(children, key=structural_record_key)
            parent_start = min(child.start for child in children)
            parent_stop = max(child.stop for child in children)

            converted.append(
                replace(
                    first,
                    feature="gene",
                    start=parent_start,
                    stop=parent_stop,
                    score=".",
                    phase=".",
                    attributes=converted_attributes(
                        first.attributes, parent_id, gene_id, locus_tag
                    ),
                )
            )
            for child in coding_children:
                converted.append(
                    replace(
                        child,
                        feature="CDS",
                        attributes=converted_attributes(
                            child.attributes,
                            f"cds{cds_count}",
                            gene_id,
                            locus_tag,
                            parent_id,
                        ),
                    )
                )
                cds_count += 1

            for child in rna_children:
                converted.append(
                    replace(
                        child,
                        phase=".",
                        attributes=converted_attributes(
                            child.attributes,
                            f"rna{rna_count}",
                            gene_id,
                            locus_tag,
                            parent_id,
                        ),
                    )
                )
                rna_count += 1

        gene_count += 1

    for unknown_count, (gene_id, record) in enumerate(unknown, start=1):
        converted.append(
            replace(
                record,
                phase="." if record.feature.lower() != "cds" else record.phase,
                attributes=converted_attributes(
                    record.attributes,
                    f"misc{unknown_count}",
                    gene_id,
                    locus_tags[gene_id],
                ),
            )
        )

    if not converted:
        raise AnnotationConversionError(
            "GTF input contains no gene, CDS, RNA, or preservable feature records"
        )
    return sorted(converted, key=output_record_key)


def write_converted(records: list[Record], output: Path) -> None:
    with output.open("w") as handle:
        handle.write("##gff-version 3\n")
        for record in records:
            handle.write(record.to_gff_line() + "\n")


def write_output(annotation: Path, output: Path) -> None:
    records, declares_gff3 = read_annotation(annotation)
    detected = annotation_format(records, declares_gff3)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    try:
        if detected == "GFF3":
            copyfile(annotation, temporary)
        else:
            write_converted(convert_gtf(records), temporary)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert GTF to GFF3, or validate and copy GFF3 input."
    )
    parser.add_argument(
        "-a", "--annotation", required=True, type=Path, help="input annotation file"
    )
    parser.add_argument(
        "-o", "--output", required=True, type=Path, help="output GFF3 annotation file"
    )
    args = parser.parse_args()

    try:
        write_output(args.annotation, args.output)
    except (AnnotationConversionError, OSError) as exc:
        parser.exit(2, f"gtf2gff3: error: {exc}\n")


if __name__ == "__main__":
    main()
