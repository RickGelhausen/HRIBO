"""
Preflight validation of HRIBO inputs.

The checks here exist to turn late, cryptic failures deep inside the workflow
(segemehl, featureCounts, bedtools, the excel writers) into a single readable
error before any job is submitted. The most common real-world cause by far is a
mismatch between the sequence identifiers used in the genome FASTA and those
used in the annotation, which nothing downstream reports intelligibly.

Only the standard library is used, so that this module can be imported from the
Snakefile at parse time, inside a plain Snakemake environment.

Author: Rick Gelhausen
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from enum import IntEnum
from pathlib import Path
from typing import Iterable, Iterator, Sequence

# Characters permitted in a nucleotide FASTA, including the IUPAC ambiguity
# codes and the gap characters that some assemblies carry.
IUPAC_NUCLEOTIDES = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn-.*")

# Feature types that HRIBO treats as protein coding. Metagene profiling, the
# ORF predictors and most of the excel output are driven entirely by these.
CDS_FEATURES = {"cds"}

# Feature types removed from the alignments by the rRNA/tRNA filtering step.
STRUCTURAL_RNA_FEATURES = {"rrna", "trna"}

MAX_REPORTED_ITEMS = 10


class Severity(IntEnum):
    """Ordered so that findings can be sorted worst-first."""

    INFO = 0
    WARNING = 1
    ERROR = 2

    @property
    def label(self) -> str:
        return self.name


@dataclass
class Finding:
    """A single validation result.

    `hint` is what the user should actually do about it, which is the part that
    was missing from the previous validation and the reason failures were hard
    to act on.
    """

    check: str
    severity: Severity
    title: str
    detail: str = ""
    hint: str = ""
    items: list[str] = field(default_factory=list)

    def format_items(self, limit: int = MAX_REPORTED_ITEMS) -> str:
        if not self.items:
            return ""
        shown = self.items[:limit]
        rendered = ", ".join(shown)
        if len(self.items) > limit:
            rendered += f", ... ({len(self.items) - limit} more)"
        return rendered

    def to_text(self) -> str:
        lines = [f"[{self.severity.label}] {self.check}: {self.title}"]
        if self.detail:
            lines.append(f"    {self.detail}")
        if self.items:
            lines.append(f"    affected: {self.format_items()}")
        if self.hint:
            lines.append(f"    hint: {self.hint}")
        return "\n".join(lines)


class ValidationReport:
    """Collects findings and decides whether the workflow may proceed."""

    def __init__(self) -> None:
        self.findings: list[Finding] = []

    def add(
        self,
        check: str,
        severity: Severity,
        title: str,
        detail: str = "",
        hint: str = "",
        items: Iterable[str] | None = None,
    ) -> None:
        self.findings.append(
            Finding(
                check=check,
                severity=severity,
                title=title,
                detail=detail,
                hint=hint,
                items=sorted(str(i) for i in items) if items else [],
            )
        )

    def error(self, check: str, title: str, **kwargs) -> None:
        self.add(check, Severity.ERROR, title, **kwargs)

    def warning(self, check: str, title: str, **kwargs) -> None:
        self.add(check, Severity.WARNING, title, **kwargs)

    def info(self, check: str, title: str, **kwargs) -> None:
        self.add(check, Severity.INFO, title, **kwargs)

    def extend(self, other: "ValidationReport") -> None:
        self.findings.extend(other.findings)

    @property
    def errors(self) -> list[Finding]:
        return [f for f in self.findings if f.severity is Severity.ERROR]

    @property
    def warnings(self) -> list[Finding]:
        return [f for f in self.findings if f.severity is Severity.WARNING]

    @property
    def ok(self) -> bool:
        return not self.errors

    def sorted_findings(self) -> list[Finding]:
        return sorted(self.findings, key=lambda f: (-f.severity, f.check))

    def to_text(self) -> str:
        if not self.findings:
            return "Input validation passed with no findings."
        body = "\n".join(f.to_text() for f in self.sorted_findings())
        summary = (
            f"{len(self.errors)} error(s), {len(self.warnings)} warning(s), "
            f"{len(self.findings)} finding(s) total."
        )
        return f"{body}\n\n{summary}"

    def raise_on_error(self, context: str = "Input validation") -> None:
        """Abort with every error at once, rather than one per run."""
        if self.ok:
            return
        blocking = "\n".join(f.to_text() for f in self.errors)
        raise ValueError(
            f"\n\n{context} failed with {len(self.errors)} error(s):\n\n"
            f"{blocking}\n\n"
            "Fix the errors above and re-run. Nothing has been executed.\n"
        )


# --------------------------------------------------------------------------
# Parsing helpers
# --------------------------------------------------------------------------


def open_maybe_gzip(path: Path):
    """Open a text file transparently, whether or not it is gzipped."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return open(path, "r", errors="replace")


@dataclass
class FastaRecord:
    identifier: str
    header: str
    length: int
    invalid_characters: set[str]


def parse_fasta(path: Path) -> list[FastaRecord]:
    """Parse a nucleotide FASTA without pulling in biopython.

    The identifier is the first whitespace-delimited token of the header, which
    is what samtools, segemehl and every GFF seqid convention use. Retaining the
    full header lets us explain mismatches caused by the description part.
    """
    records: list[FastaRecord] = []
    identifier = header = None
    length = 0
    invalid: set[str] = set()

    def flush() -> None:
        if identifier is not None:
            records.append(FastaRecord(identifier, header, length, set(invalid)))

    with open_maybe_gzip(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                identifier = header.split()[0] if header.split() else ""
                length = 0
                invalid = set()
            elif identifier is not None:
                length += len(line)
                invalid.update(set(line) - IUPAC_NUCLEOTIDES)
        flush()

    return records


@dataclass
class GffRecord:
    line_number: int
    seqid: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: str


class AnnotationParseError(Exception):
    """Raised when the annotation cannot be read as GFF/GTF at all."""


def parse_annotation(path: Path) -> tuple[list[GffRecord], bool, bool]:
    """Parse a GFF3/GTF annotation line by line.

    Returns the records, whether the file declared ``##gff-version 3`` and
    whether it carries an embedded ``##FASTA`` section. Line-based parsing is
    deliberate: reading these with ``pandas.read_csv`` breaks on the embedded
    FASTA that Prokka and NCBI routinely emit, which is one of the failure modes
    this validation exists to catch.
    """
    records: list[GffRecord] = []
    declares_gff3 = False
    has_embedded_fasta = False
    malformed: list[int] = []

    with open_maybe_gzip(path) as handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\n\r")
            if not line.strip():
                continue
            if line.startswith("##FASTA"):
                has_embedded_fasta = True
                break
            if line.startswith("#"):
                if number == 1 and line.strip() == "##gff-version 3":
                    declares_gff3 = True
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                malformed.append(number)
                continue

            seqid, source, feature, start, end, score, strand, phase, attributes = fields
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                malformed.append(number)
                continue

            records.append(
                GffRecord(
                    line_number=number,
                    seqid=seqid,
                    source=source,
                    feature=feature,
                    start=start_i,
                    end=end_i,
                    score=score,
                    strand=strand,
                    phase=phase,
                    attributes=attributes,
                )
            )

    if not records and not malformed:
        raise AnnotationParseError(
            f"No annotation records found in {path}. The file appears to be empty "
            "or to contain only comment lines."
        )

    return records, declares_gff3, has_embedded_fasta


ATTRIBUTE_GFF3_RE = re.compile(r"(?P<key>[^=;]+)=(?P<value>[^;]*)")
ATTRIBUTE_GTF_RE = re.compile(r'(?P<key>\S+)\s+"(?P<value>[^"]*)"')


def parse_attributes(attributes: str) -> dict[str, str]:
    """Parse a GFF3 or GTF attribute column into a flat dict."""
    parsed = {m.group("key").strip(): m.group("value").strip() for m in ATTRIBUTE_GTF_RE.finditer(attributes)}
    if parsed:
        return parsed
    return {m.group("key").strip(): m.group("value").strip() for m in ATTRIBUTE_GFF3_RE.finditer(attributes)}


def suggest_identifier_mapping(
    annotation_ids: Sequence[str], genome_records: Sequence[FastaRecord]
) -> list[str]:
    """Explain *why* identifiers fail to match, where the reason is recoverable.

    Mismatches are rarely arbitrary: usually the annotation uses an accession
    that appears later in the FASTA header, or differs only by a version suffix
    or by case. Naming the likely counterpart turns an opaque failure into a
    one-line fix.
    """
    suggestions: list[str] = []

    by_lower = {r.identifier.lower(): r.identifier for r in genome_records}
    by_unversioned = {r.identifier.split(".")[0]: r.identifier for r in genome_records}
    # Any whitespace-delimited token of the header, not just the first.
    by_header_token: dict[str, str] = {}
    for record in genome_records:
        for token in record.header.replace("|", " ").split():
            by_header_token.setdefault(token, record.identifier)
            by_header_token.setdefault(token.rstrip(","), record.identifier)

    for annotation_id in annotation_ids:
        candidate = None
        reason = ""
        if annotation_id.lower() in by_lower:
            candidate, reason = by_lower[annotation_id.lower()], "differs only in capitalisation"
        elif annotation_id.split(".")[0] in by_unversioned:
            candidate, reason = by_unversioned[annotation_id.split(".")[0]], "differs only in the version suffix"
        elif annotation_id in by_header_token:
            candidate, reason = (
                by_header_token[annotation_id],
                "appears in the FASTA description but not as the first token of the header",
            )
        if candidate:
            suggestions.append(f"{annotation_id!r} -> {candidate!r} ({reason})")

    return suggestions
