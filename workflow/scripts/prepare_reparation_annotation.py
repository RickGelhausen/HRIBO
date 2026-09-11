#!/usr/bin/env python
"""Build the narrow GTF dialect consumed by REPARATION 1.0.9.

REPARATION's annotation reader ignores every feature except lowercase
``transcript`` and extracts four GTF attributes with regular expressions.
HRIBO's shared annotation is deliberately GFF3, so this adapter presents each
CDS or supported structural RNA as one compatible transcript without changing
the annotation used by the rest of the workflow.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, replace
from pathlib import Path

import gff_utils


FEATURE_BIOTYPES = {
    "cds": "protein_coding",
    "ncrna": "ncRNA",
    "pseudogene": "pseudogene",
    "rrna": "rRNA",
    "srna": "sRNA",
    "trna": "tRNA",
}
EXPLICIT_TRANSCRIPTS = {"mrna", "transcript"}


class ReparationAnnotationError(ValueError):
    """The processed annotation cannot be represented safely for REPARATION."""


@dataclass(frozen=True)
class Record:
    seqid: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    attributes: dict[str, str]
    line_number: int


@dataclass(frozen=True)
class Transcript:
    seqid: str
    source: str
    start: int
    end: int
    score: str
    strand: str
    gene_id: str
    transcript_id: str
    gene_name: str
    gene_biotype: str

    def sort_key(self) -> tuple:
        return (
            self.seqid,
            self.start,
            self.end,
            self.strand,
            self.gene_id,
            self.source,
            self.score,
            self.gene_name,
            self.gene_biotype,
        )

    def to_gtf_line(self) -> str:
        attributes = (
            f'gene_id "{gtf_value(self.gene_id)}"; '
            f'transcript_id "{gtf_value(self.transcript_id)}"; '
            f'gene_name "{gtf_value(self.gene_name)}"; '
            f'gene_biotype "{gtf_value(self.gene_biotype)}";'
        )
        return "\t".join(
            (
                self.seqid,
                self.source,
                "transcript",
                str(self.start),
                str(self.end),
                self.score,
                self.strand,
                ".",
                attributes,
            )
        )


def gtf_value(value: str) -> str:
    """Protect delimiters that REPARATION's regex reader cannot escape."""

    return value.replace("\\", "%5C").replace('"', "%22")


def read_records(path: Path) -> list[Record]:
    records = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ReparationAnnotationError(
                    f"line {line_number} has {len(fields)} columns; expected nine"
                )
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ReparationAnnotationError(
                    f"line {line_number} has non-integer coordinates"
                ) from exc
            if start < 1 or end < start:
                raise ReparationAnnotationError(
                    f"line {line_number} has invalid coordinates {start}-{end}"
                )

            records.append(
                Record(
                    seqid=fields[0],
                    source=fields[1],
                    feature=fields[2],
                    start=start,
                    end=end,
                    score=fields[5],
                    strand=fields[6],
                    attributes=gff_utils.parse_attributes(fields[8]),
                    line_number=line_number,
                )
            )
    if not records:
        raise ReparationAnnotationError(f"{path} contains no annotation records")
    return records


def build_id_index(records: list[Record]) -> dict[str, Record]:
    # Some bacterial annotations reuse descriptive IDs for regulatory records
    # such as alternative TATA boxes.  Those records are not emitted for
    # REPARATION, and a duplicate cannot make parent lookup ambiguous when no
    # record refers to it.  Keep strict GFF identity checks everywhere the ID
    # can affect an emitted transcript or its ancestry.
    referenced_ids = {
        parent
        for record in records
        for parent in record.attributes.get("parent", "").split(",")
        if parent
    }
    emitted_features = set(FEATURE_BIOTYPES) | EXPLICIT_TRANSCRIPTS
    by_id = {}
    for record in records:
        identifier = record.attributes.get("id", "")
        if not identifier:
            continue
        if identifier in by_id:
            previous = by_id[identifier]
            if (
                identifier in referenced_ids
                or previous.feature.lower() in emitted_features
                or record.feature.lower() in emitted_features
            ):
                raise ReparationAnnotationError(
                    f"duplicate ID {identifier!r} on lines "
                    f"{previous.line_number} and {record.line_number}"
                )
            continue
        by_id[identifier] = record
    return by_id


def ancestors(record: Record, by_id: dict[str, Record]) -> list[Record]:
    """Return the deterministic reachable parent graph, stopping at cycles."""

    result = []
    pending = [value for value in record.attributes.get("parent", "").split(",") if value]
    seen = set()
    while pending:
        identifier = pending.pop(0)
        if identifier in seen:
            continue
        seen.add(identifier)
        parent = by_id.get(identifier)
        if parent is None:
            continue
        result.append(parent)
        pending.extend(
            value
            for value in parent.attributes.get("parent", "").split(",")
            if value and value not in seen
        )
    return result


def first_attribute(records: list[Record], *keys: str) -> str:
    for record in records:
        for key in keys:
            value = record.attributes.get(key, "")
            if value:
                return value
    return ""


def inferred_transcript_biotype(record: Record, related: list[Record]) -> str:
    feature = record.feature.lower()
    if feature == "mrna":
        return "protein_coding"
    biotype = first_attribute(
        related,
        "gene_biotype",
        "gene_type",
        "transcript_biotype",
        "transcript_type",
        "biotype",
    )
    if biotype:
        return biotype
    return "transcript"


def has_pseudogene_ancestry(records: list[Record]) -> bool:
    for record in records:
        if record.feature.lower() == "pseudogene":
            return True
        if "pseudo" in record.attributes or "pseudogene" in record.attributes:
            return True
        biotype = first_attribute(
            [record], "gene_biotype", "gene_type", "biotype"
        ).lower()
        if biotype in {"pseudo", "pseudogene"}:
            return True
    return False


def to_transcript(record: Record, by_id: dict[str, Record]) -> Transcript:
    if record.strand not in {"+", "-"}:
        raise ReparationAnnotationError(
            f"line {record.line_number} has strand {record.strand!r}; "
            "REPARATION requires '+' or '-'"
        )

    related = [record, *ancestors(record, by_id)]
    generated_id = (
        f"{record.seqid}:{record.start}-{record.end}:"
        f"{record.strand}:{record.feature.lower()}"
    )
    gene_id = (
        first_attribute(related, "gene_id", "locus_tag")
        or first_attribute(related[1:], "id")
        or first_attribute([record], "transcript_id", "id")
        or generated_id
    )
    gene_name = (
        first_attribute(related, "gene_name", "name", "gene") or gene_id
    )
    feature = record.feature.lower()
    biotype = FEATURE_BIOTYPES.get(feature)
    if feature == "cds" and has_pseudogene_ancestry(related):
        biotype = "pseudogene"
    if biotype is None:
        biotype = inferred_transcript_biotype(record, related)

    return Transcript(
        seqid=record.seqid,
        source=record.source,
        start=record.start,
        end=record.end,
        score=record.score,
        strand=record.strand,
        gene_id=gene_id,
        # REPARATION uses transcript_id only as a temporary Perl hash key. A
        # synthetic value assigned after sorting prevents duplicate source IDs
        # from silently overwriting one another inside that parser.
        transcript_id="",
        gene_name=gene_name,
        gene_biotype=biotype,
    )


def build_transcripts(records: list[Record]) -> list[Transcript]:
    by_id = build_id_index(records)
    supported = [
        record for record in records if record.feature.lower() in FEATURE_BIOTYPES
    ]

    # A pass-through GFF3 can already contain mRNA/transcript records. Keep an
    # explicit transcript only when it has no supported child: CDS/RNA geometry
    # is the unambiguous boundary REPARATION needs for ORF classification.
    supported_parent_ids = {
        parent
        for record in supported
        for parent in record.attributes.get("parent", "").split(",")
        if parent
    }
    supported = [
        record
        for record in supported
        if not (
            record.feature.lower() == "pseudogene"
            and record.attributes.get("id", "") in supported_parent_ids
        )
    ]
    for record in records:
        if record.feature.lower() not in EXPLICIT_TRANSCRIPTS:
            continue
        identifier = record.attributes.get("id", "")
        if identifier and identifier in supported_parent_ids:
            continue
        supported.append(record)

    if not supported:
        raise ReparationAnnotationError(
            "annotation has no CDS, supported RNA, pseudogene, or standalone transcript records"
        )

    source_transcripts = {}
    for record in supported:
        if record.feature.lower() != "cds":
            continue
        source_id = record.attributes.get("transcript_id", "")
        if not source_id:
            source_id = next(
                (
                    parent
                    for parent in record.attributes.get("parent", "").split(",")
                    if parent
                    and parent in by_id
                    and by_id[parent].feature.lower() in EXPLICIT_TRANSCRIPTS
                ),
                "",
            )
        if not source_id:
            continue
        previous = source_transcripts.get(source_id)
        if previous is not None:
            raise ReparationAnnotationError(
                f"CDS transcript {source_id!r} is split across lines "
                f"{previous.line_number} and {record.line_number}; "
                "REPARATION cannot represent split CDS features safely"
            )
        source_transcripts[source_id] = record

    transcripts = sorted(
        (to_transcript(record, by_id) for record in supported),
        key=Transcript.sort_key,
    )
    seen_geometry = {}
    for transcript in transcripts:
        geometry = (
            transcript.seqid,
            transcript.start,
            transcript.end,
            transcript.strand,
        )
        previous = seen_geometry.get(geometry)
        if previous is not None:
            raise ReparationAnnotationError(
                "multiple annotation records resolve to the same REPARATION "
                f"interval {transcript.seqid}:{transcript.start}-"
                f"{transcript.end}:{transcript.strand}"
            )
        seen_geometry[geometry] = transcript
    return [
        replace(
            transcript,
            transcript_id=f"hribo_reparation_transcript_{index:06d}",
        )
        for index, transcript in enumerate(transcripts, start=1)
    ]


def write_annotation(annotation: Path, output: Path) -> None:
    transcripts = build_transcripts(read_records(annotation))
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            for transcript in transcripts:
                handle.write(transcript.to_gtf_line() + "\n")
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a REPARATION-compatible transcript-only GTF."
    )
    parser.add_argument("-a", "--annotation", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    args = parser.parse_args()

    try:
        write_annotation(args.annotation, args.output)
    except (OSError, ReparationAnnotationError) as exc:
        parser.exit(2, f"prepare_reparation_annotation: error: {exc}\n")


if __name__ == "__main__":
    main()
