"""
The individual preflight checks run against HRIBO inputs.

Each check takes already-parsed inputs and returns findings; parsing lives in
lib.validation. Checks never raise for bad *data* -- they record a finding, so
that a single run reports everything wrong with the inputs instead of stopping
at the first problem.

Author: Rick Gelhausen
"""

from __future__ import annotations

import gzip
import re
import zlib
from pathlib import Path
from typing import Sequence

from lib.stages import RIBO_LIKE_METHODS, StageError, resolve_stages
from lib.validation import (
    CDS_FEATURES,
    STRUCTURAL_RNA_FEATURES,
    AnnotationParseError,
    FastaRecord,
    GffRecord,
    ValidationReport,
    parse_annotation,
    parse_attributes,
    parse_fasta,
    suggest_identifier_mapping,
)

VALID_STRANDS = {"+", "-", "."}
DNA_RE = re.compile(r"^[ACGTUNRYSWKMBDHV]+$", re.IGNORECASE)
DEEPRIBO_NUCLEOTIDES = frozenset("ACGTN")

# Ribo-seq is single-end after HRIBO's merge step; a condition needs this many
# replicates before the DESeq2-based differential expression tools are usable.
MIN_REPLICATES_FOR_DIFFEX = 2


# --------------------------------------------------------------------------
# Genome
# --------------------------------------------------------------------------


def check_genome(genome_path: Path) -> tuple[ValidationReport, list[FastaRecord]]:
    report = ValidationReport()

    if not genome_path.is_file():
        report.error(
            "GENOME_MISSING",
            "Genome file not found",
            detail=f"{genome_path} does not exist.",
            hint="Set biologySettings.genome in config/config.yaml to the path of your genome FASTA.",
        )
        return report, []

    if genome_path.stat().st_size == 0:
        report.error(
            "GENOME_EMPTY",
            "Genome file is empty",
            detail=f"{genome_path} is zero bytes.",
        )
        return report, []

    records = parse_fasta(genome_path)

    if not records:
        report.error(
            "GENOME_NOT_FASTA",
            "Genome file contains no FASTA records",
            detail=f"No '>' header line was found in {genome_path}.",
            hint="The genome must be a nucleotide FASTA file. GenBank and EMBL formats are not accepted.",
        )
        return report, records

    seen: dict[str, int] = {}
    for record in records:
        seen[record.identifier] = seen.get(record.identifier, 0) + 1
    duplicates = [name for name, count in seen.items() if count > 1]
    if duplicates:
        report.error(
            "GENOME_DUPLICATE_IDS",
            "Duplicate sequence identifiers in the genome",
            detail="samtools faidx and segemehl both refuse a FASTA with repeated identifiers.",
            hint="Rename the duplicated records so that every identifier is unique.",
            items=duplicates,
        )

    empty = [r.identifier for r in records if r.length == 0]
    if empty:
        report.error(
            "GENOME_EMPTY_SEQUENCE",
            "Genome contains records with no sequence",
            items=empty,
            hint="Remove the empty records, or replace the genome file.",
        )

    invalid = {r.identifier: sorted(r.invalid_characters) for r in records if r.invalid_characters}
    if invalid:
        report.error(
            "GENOME_INVALID_CHARACTERS",
            "Genome contains non-nucleotide characters",
            detail="Only IUPAC nucleotide codes are allowed. Amino acid FASTA files are a common mix-up here.",
            hint="Check that the file is a nucleotide and not a protein FASTA.",
            items=[f"{name} ({''.join(chars)})" for name, chars in invalid.items()],
        )

    described = [r.identifier for r in records if len(r.header.split()) > 1]
    if described:
        report.info(
            "GENOME_HEADER_DESCRIPTION",
            "Genome headers carry a description after the identifier",
            detail=(
                "Only the first whitespace-delimited token is used as the sequence identifier, "
                "which is also how the annotation seqids must be written."
            ),
            items=described,
        )

    report.info(
        "GENOME_SUMMARY",
        f"Genome has {len(records)} sequence(s), {sum(r.length for r in records):,} bp total",
        items=[f"{r.identifier} ({r.length:,} bp)" for r in records],
    )

    return report, records


def check_deepribo_genome(records: Sequence[FastaRecord]) -> ValidationReport:
    """Reject FASTA symbols that DeepRibo's sequence encoder cannot consume."""
    report = ValidationReport()
    incompatible = {
        record.identifier: sorted(
            {
                character
                for character in record.sequence_characters
                if character not in DEEPRIBO_NUCLEOTIDES
            }
        )
        for record in records
    }
    incompatible = {
        identifier: characters
        for identifier, characters in incompatible.items()
        if characters
    }
    if incompatible:
        report.error(
            "GENOME_DEEPRIBO_ALPHABET",
            "Genome contains sequence symbols unsupported by DeepRibo",
            detail=(
                "DeepRibo's case-sensitive sequence encoder accepts only uppercase "
                "A/C/G/T/N. HRIBO's other stages accept lowercase and broader IUPAC "
                "ambiguity and gap symbols, but DeepRibo would fail or lose calls "
                "while parsing these records."
            ),
            hint=(
                "Convert lowercase sequence to uppercase and replace each remaining "
                "unsupported symbol with N (or resolve it to A/C/G/T), or set "
                "predictionSettings.deepribo to 'off' if DeepRibo is not required."
            ),
            items=[
                f"{identifier} ({''.join(characters)})"
                for identifier, characters in incompatible.items()
            ],
        )
    return report


# --------------------------------------------------------------------------
# Annotation
# --------------------------------------------------------------------------


def check_annotation(annotation_path: Path) -> tuple[ValidationReport, list[GffRecord]]:
    report = ValidationReport()

    if not annotation_path.is_file():
        report.error(
            "ANNOTATION_MISSING",
            "Annotation file not found",
            detail=f"{annotation_path} does not exist.",
            hint="Set biologySettings.annotation in config/config.yaml to the path of your GFF3 or GTF file.",
        )
        return report, []

    try:
        records, declares_gff3, has_embedded_fasta, malformed = parse_annotation(
            annotation_path
        )
    except AnnotationParseError as exc:
        report.error(
            "ANNOTATION_UNPARSEABLE",
            "Annotation could not be parsed",
            detail=str(exc),
            hint="The annotation must be a tab-separated GFF3 or GTF file with nine columns.",
        )
        return report, []

    if malformed:
        report.error(
            "ANNOTATION_MALFORMED_ROWS",
            f"Annotation contains {len(malformed)} malformed row(s)",
            detail=(
                "Malformed rows cannot be represented downstream and were not included "
                "in the parsed annotation. Continuing would silently discard features."
            ),
            hint="Fix or remove every listed row; each feature must have nine tab-separated columns and integer coordinates.",
            items=malformed,
        )

    if has_embedded_fasta:
        report.error(
            "ANNOTATION_EMBEDDED_FASTA",
            "Annotation contains an embedded ##FASTA section",
            detail=(
                "Prokka and some NCBI downloads append the genome sequence to the GFF. "
                "HRIBO reads the annotation as a table and will fail on those lines."
            ),
            hint="Strip the sequence with: sed '/^##FASTA/,$d' annotation.gff > annotation.clean.gff",
        )

    bad_coordinates = [
        f"line {r.line_number} ({r.seqid}:{r.start}-{r.end})" for r in records if r.start > r.end or r.start < 1
    ]
    if bad_coordinates:
        report.error(
            "ANNOTATION_BAD_COORDINATES",
            "Annotation contains invalid coordinates",
            detail="GFF coordinates are 1-based and inclusive, so start must be >= 1 and <= end.",
            items=bad_coordinates,
        )

    bad_strands = sorted({r.strand for r in records if r.strand not in VALID_STRANDS})
    if bad_strands:
        report.error(
            "ANNOTATION_BAD_STRAND",
            "Annotation contains invalid strand values",
            detail="Column 7 must be '+', '-' or '.'.",
            items=bad_strands,
        )

    cds_records = [r for r in records if r.feature.lower() in CDS_FEATURES]
    if not cds_records:
        report.error(
            "ANNOTATION_NO_CDS",
            "Annotation contains no CDS features",
            detail=(
                "Metagene profiling, the ORF predictors and the read-count tables all select on "
                "feature type 'CDS'. Without them the workflow produces empty results."
            ),
            hint=(
                "Check column 3 of the annotation. Some Ensembl Bacteria GTF files use 'exon' plus "
                "a gene_biotype attribute rather than 'CDS'."
            ),
        )
    else:
        report.info(
            "ANNOTATION_SUMMARY",
            f"Annotation has {len(records)} feature(s), of which {len(cds_records)} are CDS",
        )

    non_triplet = [
        f"line {r.line_number} ({r.seqid}:{r.start}-{r.end}, {r.end - r.start + 1} nt)"
        for r in cds_records
        if (r.end - r.start + 1) % 3 != 0
    ]
    if non_triplet:
        report.warning(
            "ANNOTATION_CDS_NOT_TRIPLET",
            f"{len(non_triplet)} CDS feature(s) have a length that is not a multiple of three",
            detail="Reading frame and periodicity estimates are unreliable for these features.",
            items=non_triplet,
        )

    structural = [r for r in records if r.feature.lower() in STRUCTURAL_RNA_FEATURES]
    if not structural:
        report.warning(
            "ANNOTATION_NO_RRNA",
            "Annotation contains no rRNA or tRNA features",
            detail=(
                "The rRNA/tRNA filtering step builds its BED file by selecting those feature types. "
                "With none present the filter is a no-op and structural RNA reads stay in the data, "
                "which typically dominate a Ribo-seq library."
            ),
            hint="Verify that your annotation includes rRNA and tRNA features before trusting the read counts.",
        )
    else:
        report.info(
            "ANNOTATION_RRNA_SUMMARY",
            f"Annotation has {len(structural)} rRNA/tRNA feature(s) available for filtering",
        )

    _check_annotation_attributes(report, cds_records, declares_gff3)

    return report, records


def _check_annotation_attributes(
    report: ValidationReport, cds_records: Sequence[GffRecord], declares_gff3: bool
) -> None:
    """Attribute-level checks for the fields HRIBO's output tables depend on."""
    if not cds_records:
        return

    identifiers: dict[str, list[int]] = {}
    missing_locus_tag = []

    for record in cds_records:
        attributes = parse_attributes(record.attributes)
        identifier = attributes.get("ID") or attributes.get("gene_id")
        if identifier:
            identifiers.setdefault(identifier, []).append(record.line_number)
        if not (attributes.get("locus_tag") or attributes.get("gene_id")):
            missing_locus_tag.append(record.line_number)

    duplicates = {k: v for k, v in identifiers.items() if len(v) > 1}
    if duplicates:
        report.warning(
            "ANNOTATION_DUPLICATE_IDS",
            f"{len(duplicates)} feature identifier(s) are used more than once",
            detail="Read counts are keyed by feature identifier, so duplicates are merged or overwritten.",
            items=[f"{k} (lines {', '.join(str(n) for n in v)})" for k, v in duplicates.items()],
        )

    if missing_locus_tag:
        report.warning(
            "ANNOTATION_NO_LOCUS_TAG",
            f"{len(missing_locus_tag)} CDS feature(s) have neither locus_tag nor gene_id",
            detail="These features appear unnamed in the excel output tables.",
            items=[f"line {n}" for n in missing_locus_tag],
        )

    if not identifiers:
        report.warning(
            "ANNOTATION_NO_IDENTIFIERS",
            "No CDS feature carries an ID or gene_id attribute",
            detail=(
                f"The annotation looks like {'GFF3' if declares_gff3 else 'GTF'}. "
                "Attribute parsing may be picking up the wrong format."
            ),
            hint="Check that column 9 uses 'key=value;' (GFF3) or 'key \"value\";' (GTF) consistently.",
        )


# --------------------------------------------------------------------------
# Genome / annotation consistency -- the check this module mainly exists for
# --------------------------------------------------------------------------


def check_reference_consistency(
    genome_records: Sequence[FastaRecord], gff_records: Sequence[GffRecord]
) -> ValidationReport:
    report = ValidationReport()

    if not genome_records or not gff_records:
        return report

    genome_ids = {r.identifier for r in genome_records}
    genome_lengths = {r.identifier: r.length for r in genome_records}
    annotation_ids = {r.seqid for r in gff_records}

    orphaned = sorted(annotation_ids - genome_ids)
    if orphaned:
        suggestions = suggest_identifier_mapping(orphaned, genome_records)
        hint = (
            "Rename the sequences so that both files agree. "
            "The identifier is the first whitespace-delimited token of the FASTA header."
        )
        if suggestions:
            hint = "Likely intended matches:\n      " + "\n      ".join(suggestions) + "\n    " + hint
        report.error(
            "REFERENCE_ID_MISMATCH",
            f"{len(orphaned)} annotation sequence identifier(s) are absent from the genome",
            detail=(
                "Every seqid in column 1 of the annotation must name a sequence in the genome FASTA. "
                "Where they disagree, mapping and read counting silently produce empty results or "
                "fail much later with an unrelated-looking error.\n"
                f"    genome has: {', '.join(sorted(genome_ids))}\n"
                f"    annotation has: {', '.join(sorted(annotation_ids))}"
            ),
            hint=hint,
            items=orphaned,
        )

    unannotated = sorted(genome_ids - annotation_ids)
    if unannotated:
        report.warning(
            "REFERENCE_UNANNOTATED_SEQUENCE",
            f"{len(unannotated)} genome sequence(s) have no annotated features",
            detail="Reads mapping to these sequences are counted as mapped but never assigned to a feature.",
            items=unannotated,
        )

    out_of_bounds = [
        f"line {r.line_number} ({r.seqid}:{r.start}-{r.end}, sequence is {genome_lengths[r.seqid]:,} bp)"
        for r in gff_records
        if r.seqid in genome_lengths and r.end > genome_lengths[r.seqid]
    ]
    if out_of_bounds:
        report.error(
            "REFERENCE_COORDINATE_OUT_OF_BOUNDS",
            f"{len(out_of_bounds)} annotated feature(s) extend past the end of their sequence",
            detail="This usually means the annotation was built against a different assembly version.",
            hint="Confirm that the genome and the annotation come from the same assembly release.",
            items=out_of_bounds,
        )

    return report


# --------------------------------------------------------------------------
# Sample sheet
# --------------------------------------------------------------------------


def check_sample_sheet(samples) -> ValidationReport:
    """Cross-row checks that a JSON schema cannot express."""
    report = ValidationReport()

    if len(samples) == 0:
        report.error("SAMPLES_EMPTY", "Sample sheet contains no rows")
        return report

    keys = samples[["method", "condition", "replicate"]].astype(str).agg("-".join, axis=1)
    duplicated = sorted(set(keys[keys.duplicated()]))
    if duplicated:
        report.error(
            "SAMPLES_DUPLICATE_KEY",
            "Sample sheet contains duplicate method-condition-replicate combinations",
            detail=(
                "Every intermediate file is named <method>-<condition>-<replicate>, so duplicate "
                "keys would have different libraries overwrite each other."
            ),
            items=duplicated,
        )

    for column in ("method", "condition"):
        offending = sorted(
            {value for value in samples[column].astype(str) if not value.isalnum()}
        )
        if offending:
            report.error(
                f"SAMPLES_{column.upper()}_NOT_ALNUM",
                f"Column '{column}' contains non-alphanumeric values",
                detail=(
                    "These values become part of both filenames and Snakemake wildcards. A '-' in "
                    "particular makes <method>-<condition>-<replicate> ambiguous and breaks "
                    "contrast parsing."
                ),
                items=offending,
            )

    fastq_columns = [c for c in ("fastqFile", "fastqFile2") if c in samples.columns]
    used: dict[str, list[str]] = {}
    for _, row in samples.iterrows():
        label = f"{row['method']}-{row['condition']}-{row['replicate']}"
        for column in fastq_columns:
            value = row.get(column)
            if isinstance(value, str) and value.strip():
                used.setdefault(value.strip(), []).append(label)
    reused = {path: labels for path, labels in used.items() if len(labels) > 1}
    if reused:
        report.warning(
            "SAMPLES_REUSED_FASTQ",
            f"{len(reused)} fastq file(s) are referenced by more than one library",
            detail="Usually a copy-paste error in the sample sheet.",
            items=[f"{path} (used by {', '.join(labels)})" for path, labels in reused.items()],
        )

    return report


def _fastq_format_problem(path: Path) -> str | None:
    """Stream ``path`` and return its first FASTQ structure error, if any."""
    first_problem: str | None = None
    record_number = 0

    # This deliberately scans the complete decompressed file, so preflight I/O
    # is linear in FASTQ size. Reaching EOF is the only way gzip can verify its
    # trailer/CRC, and later records need the same validation as the first one;
    # only the current four-line record is retained in memory.
    with gzip.open(path, "rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                break

            record_number += 1
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()

            missing = [
                name
                for name, line in (
                    ("sequence line", sequence),
                    ("separator line", separator),
                    ("quality line", quality),
                )
                if not line
            ]
            if missing:
                if first_problem is None:
                    first_problem = (
                        f"record {record_number} is incomplete (missing {', '.join(missing)})"
                    )
                break

            if first_problem is not None:
                continue
            if not header.startswith(b"@"):
                first_problem = f"record {record_number} header line must start with '@'"
            elif not separator.startswith(b"+"):
                first_problem = f"record {record_number} separator line must start with '+'"
            else:
                sequence_length = len(sequence.rstrip(b"\r\n"))
                quality_length = len(quality.rstrip(b"\r\n"))
                if sequence_length != quality_length:
                    first_problem = (
                        f"record {record_number} sequence/quality length mismatch "
                        f"({sequence_length} != {quality_length})"
                    )

    if record_number == 0:
        return "contains no FASTQ records"
    return first_problem


def check_fastq_files(samples) -> ValidationReport:
    """Verify referenced FASTQs, including complete gzip and record integrity."""
    report = ValidationReport()

    missing: list[str] = []
    unreadable: list[str] = []
    empty: list[str] = []
    not_gzip: list[str] = []
    malformed: list[str] = []
    validated_paths: set[Path] = set()

    for _, row in samples.iterrows():
        for column in ("fastqFile", "fastqFile2"):
            value = row.get(column)
            if not isinstance(value, str) or not value.strip():
                continue
            path = Path(value.strip())
            if not path.is_file():
                missing.append(f"{value} ({column} of {row['method']}-{row['condition']}-{row['replicate']})")
                continue
            if path.stat().st_size == 0:
                empty.append(str(path))
                continue
            # Reuse is reported by check_sample_sheet; avoid decompressing the
            # same exact multi-gigabyte path once for every referencing row.
            if path in validated_paths:
                continue
            validated_paths.add(path)
            try:
                with path.open("rb") as handle:
                    if handle.read(2) != b"\x1f\x8b":
                        not_gzip.append(str(path))
                        continue
                problem = _fastq_format_problem(path)
                if problem is not None:
                    malformed.append(f"{path} ({problem})")
            except (OSError, EOFError, zlib.error) as exc:
                reason = f"{exc.__class__.__name__}: {exc}" if str(exc) else exc.__class__.__name__
                unreadable.append(f"{path} ({reason})")

    if missing:
        report.error(
            "FASTQ_MISSING",
            f"{len(missing)} fastq file(s) referenced by the sample sheet do not exist",
            hint="Paths are resolved relative to the working directory passed via --directory.",
            items=missing,
        )
    if empty:
        report.error("FASTQ_EMPTY", f"{len(empty)} fastq file(s) are empty", items=empty)
    if unreadable:
        report.error(
            "FASTQ_UNREADABLE",
            f"{len(unreadable)} fastq file(s) could not be read",
            detail=(
                "For gzip inputs, truncation or a CRC error usually means a failed download "
                "or an interrupted copy."
            ),
            items=unreadable,
        )
    if not_gzip:
        report.error(
            "FASTQ_NOT_GZIP",
            f"{len(not_gzip)} fastq file(s) are not gzip-compressed",
            detail=(
                "HRIBO stages every input under a '.fastq.gz' name, so the file contents "
                "must be gzip-compressed regardless of the source filename extension."
            ),
            items=not_gzip,
        )
    if malformed:
        report.error(
            "FASTQ_MALFORMED",
            f"{len(malformed)} file(s) do not contain valid four-line FASTQ records",
            detail=(
                "Each record needs an '@' header, a sequence, a '+' separator, and an "
                "equal-length quality string."
            ),
            items=malformed,
        )

    return report


# --------------------------------------------------------------------------
# Config semantics
# --------------------------------------------------------------------------


def check_config_semantics(config: dict, samples) -> ValidationReport:
    """Checks that depend on the config and the sample sheet together.

    Structural validation of the config is handled by the JSON schema; what is
    left here is everything the schema cannot see, such as whether the requested
    contrasts actually exist among the sampled conditions.
    """
    report = ValidationReport()

    biology = config.get("biologySettings", {})
    diffex = config.get("differentialExpressionSettings", {})
    metagene = config.get("metageneSettings", {})

    conditions = set(samples["condition"].astype(str))
    methods = set(samples["method"].astype(str))

    _check_adapters(report, biology)

    stages = _check_stages(report, config, methods)

    if "differential_expression" in stages:
        _check_diffex_feasibility(report, diffex, samples, conditions, methods)

    if "pca" in stages and len(samples.index) < 2:
        report.error(
            "PCA_TOO_FEW_LIBRARIES",
            "PCA requires at least two libraries",
            detail="A single observation has no between-library variance to decompose.",
            hint="Add another library or drop the 'pca' stage.",
            items=[
                f"{row.method}-{row.condition}-{row.replicate}"
                for row in samples.itertuples(index=False)
            ],
        )

    if "metagene" in stages or "tis_advisor" in stages:
        _check_metagene_settings(report, metagene, methods)

    return report


def _check_stages(report: ValidationReport, config: dict, methods: set[str]) -> list[str]:
    """Resolve the requested stages, reporting what cannot be run.

    Returns the stages that will actually be built, so the checks below only
    complain about settings this run depends on.
    """
    has_ribo = "RIBO" in methods
    has_ribo_like = bool(RIBO_LIKE_METHODS & methods)
    try:
        requested = resolve_stages(config, has_ribo=True)
    except StageError as exc:
        report.error(
            "CONFIG_UNKNOWN_STAGE",
            "The requested workflow stages could not be interpreted",
            detail=str(exc),
        )
        return []

    stages = resolve_stages(
        config,
        has_ribo=has_ribo,
        has_ribo_like=has_ribo_like,
    )

    skipped = [stage for stage in requested if stage not in stages]
    if skipped:
        report.warning(
            "STAGES_NEED_RIBO",
            "Requested stages were skipped because their required ribosome-profiling library type is absent",
            detail="They are skipped; the rest of the requested stages still run.",
            items=skipped,
        )

    return stages


def _check_adapters(report: ValidationReport, biology: dict) -> None:
    invalid = []
    for key in ("adapterS3", "adapterS5", "adapterP3R1", "adapterP5R1", "adapterP3R2", "adapterP5R2"):
        value = biology.get(key, "")
        if not isinstance(value, str) or not value.strip():
            continue
        for adapter in value.split(","):
            adapter = adapter.strip()
            if adapter and not DNA_RE.match(adapter):
                invalid.append(f"{key}={adapter!r}")

    if invalid:
        report.error(
            "CONFIG_INVALID_ADAPTER",
            "Adapter sequences contain non-nucleotide characters",
            detail="cutadapt rejects adapters that are not nucleotide sequences.",
            items=invalid,
        )


def _check_diffex_feasibility(
    report: ValidationReport, diffex: dict, samples, conditions: set[str], methods: set[str]
) -> None:
    contrasts = diffex.get("contrasts") or []
    ribo_conditions = set(samples.loc[samples["method"] == "RIBO", "condition"])
    rna_conditions = set(samples.loc[samples["method"] == "RNA", "condition"])
    matched_conditions = ribo_conditions & rna_conditions

    if not contrasts and len(matched_conditions) < 2:
        report.error(
            "DIFFEX_SINGLE_CONDITION",
            "Differential expression needs at least two conditions with both RIBO and RNA libraries",
            detail=(
                "Automatic contrasts are derived only from matched RIBO/RNA conditions; "
                "conditions belonging only to another assay are not comparable here."
            ),
            hint="Add a second matched condition, configure an explicit valid contrast, or drop the 'differential_expression' stage.",
            items=sorted(matched_conditions),
        )

    unknown = []
    self_contrasts = []
    selected_conditions = set() if contrasts else set(matched_conditions)
    for contrast in contrasts:
        parts = str(contrast).split("-")
        if len(parts) != 2:
            report.error(
                "DIFFEX_MALFORMED_CONTRAST",
                f"Contrast {contrast!r} is not of the form 'treated-untreated'",
                detail="Condition names may not contain '-', because it separates the two sides of a contrast.",
            )
            continue
        if parts[0] == parts[1]:
            self_contrasts.append(str(contrast))
        unknown.extend(p for p in parts if p not in conditions)
        selected_conditions.update(p for p in parts if p in conditions)
    if self_contrasts:
        report.error(
            "DIFFEX_SELF_CONTRAST",
            "Differential-expression contrasts must compare two different conditions",
            detail=(
                "A self-contrast duplicates the same count columns on both sides and "
                "produces a one-level statistical design rather than a comparison."
            ),
            hint="Remove the self-contrast or replace one side with a different condition.",
            items=sorted(set(self_contrasts)),
        )
    if unknown:
        report.error(
            "DIFFEX_UNKNOWN_CONDITION",
            "Configured contrasts reference conditions that are not in the sample sheet",
            hint=f"Available conditions: {', '.join(sorted(conditions))}",
            items=sorted(set(unknown)),
        )

    # xtail, riborex and deltaTE all compute translational efficiency, which
    # needs a matched RNA-seq library for every Ribo-seq library.
    if "RIBO" in methods and "RNA" not in methods:
        report.error(
            "DIFFEX_NO_RNA",
            "Differential expression is enabled but the sample sheet contains no RNA libraries",
            detail=(
                "xtail, riborex and deltaTE all compare Ribo-seq against RNA-seq to derive "
                "translational efficiency; none of them can run on Ribo-seq alone."
            ),
            hint="Add the matching RNA-seq libraries, or turn differential expression off.",
        )

    for method in ("RIBO", "RNA"):
        if method not in methods:
            continue
        for condition in sorted(selected_conditions):
            count = len(samples[(samples["method"] == method) & (samples["condition"] == condition)])
            if count == 0:
                report.error(
                    "DIFFEX_MISSING_GROUP",
                    f"No {method} library for condition {condition!r}",
                    detail="Every condition entering a contrast needs both a RIBO and an RNA library.",
                )
            elif count < MIN_REPLICATES_FOR_DIFFEX:
                report.error(
                    "DIFFEX_SINGLE_REPLICATE",
                    f"Only {count} {method} replicate for condition {condition!r}",
                    detail=(
                        "DESeq2, which underlies deltaTE and riborex, cannot estimate "
                        "dispersion from a single replicate. HRIBO does not schedule "
                        "differential-expression tools for an unsupported design."
                    ),
                    hint="Provide at least two biological replicates for every RIBO and RNA condition.",
                )

    if {"RIBO", "RNA"} <= methods:
        for condition in sorted(selected_conditions):
            ribo_count = len(
                samples[
                    (samples["method"] == "RIBO")
                    & (samples["condition"] == condition)
                ]
            )
            rna_count = len(
                samples[
                    (samples["method"] == "RNA")
                    & (samples["condition"] == condition)
                ]
            )
            if ribo_count and rna_count and ribo_count != rna_count:
                report.error(
                    "DIFFEX_UNMATCHED_REPLICATES",
                    f"Condition {condition!r} has {ribo_count} RIBO but {rna_count} RNA replicates",
                    detail=(
                        "Riborex reuses the RIBO condition vector for the RNA matrix, so "
                        "unequal table widths would mislabel samples or fail inside R."
                    ),
                    hint="Provide the same number of RIBO and RNA biological replicates for each contrasted condition.",
                )


def _check_metagene_settings(report: ValidationReport, metagene: dict, methods: set[str]) -> None:
    if not metagene:
        return

    ribo_like = {"RIBO", "TIS", "TTS"} & methods
    if not ribo_like:
        report.warning(
            "METAGENE_NO_RIBO",
            "No RIBO, TIS or TTS libraries present, so no metagene profiling will be produced",
        )

    positions_in_orf = metagene.get("positionsInORF")
    length_cutoff = metagene.get("lengthCutoff")
    if isinstance(positions_in_orf, int) and isinstance(length_cutoff, int):
        if length_cutoff < positions_in_orf:
            report.info(
                "METAGENE_LENGTH_CUTOFF",
                f"lengthCutoff ({length_cutoff}) is below positionsInORF ({positions_in_orf})",
                detail=(
                    "Genes shorter than positionsInORF are dropped regardless, so the effective "
                    f"minimum gene length is {positions_in_orf} nt."
                ),
            )


def check_metagene_annotation_coverage(
    config: dict,
    gff_records: Sequence[GffRecord],
    methods: set[str] | None = None,
) -> ValidationReport:
    """Warn when the metagene window filters away almost the whole annotation."""
    report = ValidationReport()

    try:
        if methods is None:
            stages = resolve_stages(config)
        else:
            stages = resolve_stages(
                config,
                has_ribo="RIBO" in methods,
                has_ribo_like=bool(RIBO_LIKE_METHODS & methods),
            )
    except StageError:
        return report  # reported by check_config_semantics
    if "metagene" not in stages and "tis_advisor" not in stages:
        return report

    metagene = config.get("metageneSettings", {})
    positions_in_orf = metagene.get("positionsInORF")
    if not isinstance(positions_in_orf, int) or not gff_records:
        return report

    cds = [r for r in gff_records if r.feature.lower() in CDS_FEATURES]
    if not cds:
        return report

    surviving = [r for r in cds if (r.end - r.start + 1) >= positions_in_orf]
    fraction = len(surviving) / len(cds)

    if not surviving:
        report.error(
            "METAGENE_NO_GENES_SURVIVE",
            f"No CDS is at least positionsInORF ({positions_in_orf} nt) long",
            detail="Metagene profiling would run on an empty gene set and produce empty plots.",
            hint="Lower metageneSettings.positionsInORF.",
        )
    elif fraction < 0.1:
        report.warning(
            "METAGENE_FEW_GENES_SURVIVE",
            f"Only {len(surviving)} of {len(cds)} CDS ({fraction:.1%}) are long enough for the metagene window",
            detail=f"positionsInORF is {positions_in_orf} nt.",
            hint="Consider lowering metageneSettings.positionsInORF to retain more genes.",
        )
    else:
        report.info(
            "METAGENE_GENES_SURVIVE",
            f"{len(surviving)} of {len(cds)} CDS ({fraction:.1%}) pass the metagene length filter",
        )

    return report
