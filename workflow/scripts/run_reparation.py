#!/usr/bin/env python3
"""Run REPARATION behind a validated, transactional output boundary."""

import argparse
from contextlib import contextmanager
import csv
import fcntl
import hashlib
import math
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path


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
OFFSET_HEADER = ["length", "p_offset"]
PDF_NAMES = (
    "metagene_profile.pdf",
    "PR_and_ROC_curve.pdf",
    "variable_importance.pdf",
    "S_Curve.pdf",
)
SAFE_PATH = re.compile(r"^[A-Za-z0-9_./-]+$")
LOCUS = re.compile(r"^(.+):(\d+)-(\d+)$")
INTEGER = re.compile(r"^[+-]?\d+$")
CODON = re.compile(r"^[ACGT]{3}$", re.IGNORECASE)
FASTA_SEQUENCE = re.compile(r"^[A-Za-z*.-]+$")
GENOME_SEQUENCE = re.compile(r"^[ACGTRYSWKMBDHVNU]+$", re.IGNORECASE)
FASTA_HEADER = re.compile(
    r"^generic\|(.+)\|start codon:([ACGT]{3}) strand:([+-]) length:(\d+)$",
    re.IGNORECASE,
)
COMPLETION_RECEIPT = ".complete"
COMPLETION_RECEIPT_TEMPORARY = ".complete.tmp"
COMPLETION_RECEIPT_CONTENT = b"HRIBO REPARATION publication v1\n"
STAGING_LOCATOR_VERSION = "HRIBO REPARATION staging locator v1"
STAGING_MARKER = ".hribo_reparation_stage"
STAGING_MARKER_VERSION = "HRIBO REPARATION staging v1"


class ArtifactError(RuntimeError):
    """A REPARATION artifact is absent or violates its file contract."""


def _arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", default="reparation.pl")
    parser.add_argument("--genome", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bai", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", required=True, type=int)
    return parser.parse_args()


def _input_file(raw_path, label):
    path = Path(raw_path).expanduser().resolve()
    if not path.is_file():
        raise ArtifactError("{} input is not a file: {}".format(label, raw_path))
    return path


def _output_directory(raw_path):
    path = Path(os.path.abspath(os.path.expanduser(raw_path)))
    if path.parent == path:
        raise ArtifactError("refusing to use the filesystem root as --output-dir")
    expected_parent = Path.cwd() / "reparation"
    if path.parent != expected_parent:
        raise ArtifactError(
            "--output-dir must be a direct child of the workflow's reparation "
            "directory: {}".format(expected_parent)
        )
    if os.path.lexists(str(expected_parent)) and expected_parent.is_symlink():
        raise ArtifactError(
            "refusing to publish through a symlinked reparation directory: {}".format(
                expected_parent
            )
        )
    if os.path.lexists(str(path)) and (path.is_symlink() or not path.is_dir()):
        raise ArtifactError(
            "existing --output-dir must be a real directory, not a file or symlink: {}".format(
                raw_path
            )
        )
    return path


def _validate_output_inputs(output_dir, inputs):
    for label, source in inputs.items():
        if output_dir == source or output_dir in source.parents:
            raise ArtifactError(
                "refusing to replace --output-dir because it contains the {} "
                "input: {}".format(label, source)
            )


def _safe_temporary_base():
    # tempfile.gettempdir() honors a scheduler-provided TMPDIR. Fall back to
    # /tmp only when that location cannot safely be passed to the legacy tool.
    candidates = [Path(tempfile.gettempdir()), Path("/tmp")]
    for candidate in candidates:
        try:
            resolved = candidate.resolve()
        except OSError:
            continue
        if (
            resolved.is_dir()
            and os.access(str(resolved), os.W_OK | os.X_OK)
            and SAFE_PATH.fullmatch(str(resolved))
        ):
            mode = resolved.stat().st_mode
            # A deterministic child can only be protected from replacement by
            # other users when the temporary base is private or sticky. This
            # includes scheduler-owned 0700 scratch directories and /tmp.
            if not mode & 0o022 or mode & stat.S_ISVTX:
                return str(resolved)
    raise ArtifactError("no shell-safe temporary directory is available")


def _staging_digest(output_dir):
    return hashlib.sha256(os.fsencode(str(output_dir))).hexdigest()


def _staging_namespace(temporary_base):
    return Path(temporary_base) / "hribo_reparation_{}".format(os.getuid())


def _staging_path(output_dir, temporary_base):
    return _staging_namespace(temporary_base) / "output_{}".format(
        _staging_digest(output_dir)
    )


def _staging_marker_content(output_dir):
    return "{}\n{}\n".format(
        STAGING_MARKER_VERSION, _staging_digest(output_dir)
    ).encode("ascii")


def _stage_inputs(stage_root, genome, gtf, database, bam, bai):
    inputs = stage_root / "inputs"
    inputs.mkdir()

    aliases = {
        "genome": inputs / "genome.fa",
        "gtf": inputs / "annotation.gtf",
        "database": inputs / "protein_db.fasta",
        "bam": inputs / "reads.bam",
        "bai": inputs / "reads.bam.bai",
    }
    aliases["genome"].symlink_to(genome)
    aliases["gtf"].symlink_to(gtf)
    aliases["bam"].symlink_to(bam)
    aliases["bai"].symlink_to(bai)
    # makeblastdb names its sidecars from the alias passed with -in. A symlink
    # therefore avoids copying Swiss-Prot while keeping all mutable indices in
    # this per-run staging directory.
    aliases["database"].symlink_to(database)
    return aliases


def _require_file(path, label, allow_empty=False):
    if path.is_symlink() or not path.is_file():
        raise ArtifactError("REPARATION did not produce required {}: {}".format(label, path.name))
    if not allow_empty and path.stat().st_size == 0:
        raise ArtifactError("REPARATION produced an empty {}: {}".format(label, path.name))


def _finite_number(value, field, line_number, minimum=None, maximum=None):
    try:
        number = float(value)
    except ValueError:
        raise ArtifactError(
            "Predicted_ORFs.txt line {} has a non-numeric {}".format(line_number, field)
        )
    if not math.isfinite(number):
        raise ArtifactError(
            "Predicted_ORFs.txt line {} has a non-finite {}".format(line_number, field)
        )
    if minimum is not None and number < minimum:
        raise ArtifactError(
            "Predicted_ORFs.txt line {} has {} below {}".format(
                line_number, field, minimum
            )
        )
    if maximum is not None and number > maximum:
        raise ArtifactError(
            "Predicted_ORFs.txt line {} has {} above {}".format(
                line_number, field, maximum
            )
        )
    return number


def _integer_or_na(value, field, line_number):
    if value == "NA":
        return
    if not INTEGER.fullmatch(value):
        raise ArtifactError(
            "Predicted_ORFs.txt line {} has an invalid {}".format(line_number, field)
        )


def _reference_lengths(path):
    """Return FASTA record lengths without loading the genome into memory."""
    lengths = {}
    identifier = None
    sequence_length = 0
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    if sequence_length == 0:
                        raise ArtifactError(
                            "genome FASTA record {!r} has no sequence".format(identifier)
                        )
                    lengths[identifier] = sequence_length
                fields = line[1:].split()
                if not fields:
                    raise ArtifactError(
                        "genome FASTA line {} has an empty header".format(line_number)
                    )
                identifier = fields[0]
                if identifier in lengths:
                    raise ArtifactError(
                        "genome FASTA has duplicate record {!r}".format(identifier)
                    )
                sequence_length = 0
            else:
                if identifier is None:
                    raise ArtifactError(
                        "genome FASTA has sequence before its first header"
                    )
                if not GENOME_SEQUENCE.fullmatch(line):
                    raise ArtifactError(
                        "genome FASTA line {} contains non-IUPAC sequence data".format(
                            line_number
                        )
                    )
                sequence_length += len(line)
    if identifier is None:
        raise ArtifactError("genome FASTA has no records")
    if sequence_length == 0:
        raise ArtifactError(
            "genome FASTA record {!r} has no sequence".format(identifier)
        )
    lengths[identifier] = sequence_length
    return lengths


def _validate_orf_table(path, reference_lengths):
    _require_file(path, "annotated ORF table")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != ORF_HEADER:
        raise ArtifactError("Predicted_ORFs.txt does not have the exact 13-column header")

    records = []
    identities = set()
    for line_number, row in enumerate(rows[1:], 2):
        if len(row) != len(ORF_HEADER):
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has {} columns; expected 13".format(
                    line_number, len(row)
                )
            )
        locus_match = LOCUS.fullmatch(row[0])
        if not locus_match:
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has an invalid ORF_locus".format(line_number)
            )
        start = int(locus_match.group(2))
        stop = int(locus_match.group(3))
        contig = locus_match.group(1)
        if start < 1 or stop < start:
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has an invalid coordinate range".format(
                    line_number
                )
            )
        if row[1] not in ("+", "-"):
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has an invalid strand".format(line_number)
            )
        if contig not in reference_lengths:
            raise ArtifactError(
                "Predicted_ORFs.txt line {} refers to unknown genome contig {!r}".format(
                    line_number, contig
                )
            )
        complete_start = start if row[1] == "+" else start - 3
        complete_stop = stop + 3 if row[1] == "+" else stop
        if complete_start < 1 or complete_stop > reference_lengths[contig]:
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has stop-inclusive interval {}-{} "
                "outside contig {!r} (length {})".format(
                    line_number,
                    complete_start,
                    complete_stop,
                    contig,
                    reference_lengths[contig],
                )
            )
        if (
            not INTEGER.fullmatch(row[2])
            or int(row[2]) != stop - start + 1
            or int(row[2]) % 3 != 0
        ):
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has a non-coding length or a length "
                "inconsistent with its locus".format(line_number)
            )
        if not CODON.fullmatch(row[3]):
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has an invalid start codon".format(line_number)
            )
        _finite_number(row[4], "ribo_count", line_number, minimum=0)
        _finite_number(row[5], "ribo_rpkm", line_number, minimum=0)
        _finite_number(row[6], "ribo_coverage", line_number, minimum=0, maximum=1)
        _finite_number(row[7], "SD_score", line_number)
        _integer_or_na(row[8], "SD_pos", line_number)
        _finite_number(row[9], "prob", line_number, minimum=0, maximum=1)
        if not row[10]:
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has an empty ORF_type".format(line_number)
            )
        if not row[11]:
            raise ArtifactError(
                "Predicted_ORFs.txt line {} has an empty Reference".format(line_number)
            )
        _integer_or_na(row[12], "Distance_from_aTIS", line_number)
        identity = (row[0], row[1])
        if identity in identities:
            raise ArtifactError(
                "Predicted_ORFs.txt has duplicate ORF {} on strand {}".format(
                    row[0], row[1]
                )
            )
        identities.add(identity)
        records.append(
            {
                "locus": row[0],
                "contig": contig,
                "start": start,
                "stop": stop,
                "strand": row[1],
                "length": int(row[2]),
                "codon": row[3].upper(),
            }
        )
    return records


def _validate_bed(path, expected_records):
    _require_file(path, "BED file")
    with path.open(encoding="utf-8") as handle:
        lines = [line.rstrip("\r\n") for line in handle]
    if not lines or not lines[0].startswith("track type=bed "):
        raise ArtifactError("Predicted_ORFs.bed does not start with a BED track header")
    records = lines[1:]
    if len(records) != len(expected_records):
        raise ArtifactError(
            "Predicted_ORFs.bed has {} records; expected {}".format(
                len(records), len(expected_records)
            )
        )
    expected = {
        (record["locus"], record["strand"]): record
        for record in expected_records
    }
    seen = set()
    for line_number, line in enumerate(records, 2):
        fields = line.split("\t")
        if len(fields) != 9:
            raise ArtifactError(
                "Predicted_ORFs.bed line {} has {} columns; expected 9".format(
                    line_number, len(fields)
                )
            )
        try:
            start = int(fields[1])
            end = int(fields[2])
            score = int(fields[4])
            thick_start = int(fields[6])
            thick_end = int(fields[7])
        except ValueError:
            raise ArtifactError(
                "Predicted_ORFs.bed line {} has non-integer coordinates or score".format(
                    line_number
                )
            )
        if not fields[0] or start < 0 or end <= start:
            raise ArtifactError(
                "Predicted_ORFs.bed line {} has an invalid interval".format(line_number)
            )
        if not 0 <= score <= 1000 or fields[5] not in ("+", "-"):
            raise ArtifactError(
                "Predicted_ORFs.bed line {} has an invalid score or strand".format(
                    line_number
                )
            )
        if not start <= thick_start <= thick_end <= end:
            raise ArtifactError(
                "Predicted_ORFs.bed line {} has invalid thick coordinates".format(
                    line_number
                )
            )
        identity = (fields[3], fields[5])
        if identity not in expected:
            raise ArtifactError(
                "Predicted_ORFs.bed line {} does not match an ORF table locus and "
                "strand".format(line_number)
            )
        if identity in seen:
            raise ArtifactError(
                "Predicted_ORFs.bed contains duplicate record {} on strand {}".format(
                    fields[3], fields[5]
                )
            )
        seen.add(identity)
        record = expected[identity]
        if record["strand"] == "+":
            expected_coordinates = (
                record["start"] - 1,
                record["stop"] + 3,
                record["start"] - 1,
                record["stop"],
            )
        else:
            expected_coordinates = (
                record["start"] - 4,
                record["stop"],
                record["start"] - 1,
                record["stop"],
            )
        coordinates = (start, end, thick_start, thick_end)
        if fields[0] != record["contig"] or coordinates != expected_coordinates:
            raise ArtifactError(
                "Predicted_ORFs.bed line {} has a contig or coordinates "
                "inconsistent with the ORF table".format(line_number)
            )
    if seen != set(expected):
        raise ArtifactError("Predicted_ORFs.bed records do not match the ORF table")


def _validate_fasta(path, expected_records):
    _require_file(path, "protein FASTA", allow_empty=not expected_records)
    fasta_records = []
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(line) == 1:
                    raise ArtifactError(
                        "Predicted_ORFs.fasta line {} has an empty header".format(
                            line_number
                        )
                    )
                if header is not None and not sequence:
                    raise ArtifactError("Predicted_ORFs.fasta contains an empty record")
                if header is not None:
                    fasta_records.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            else:
                if header is None or not FASTA_SEQUENCE.fullmatch(line):
                    raise ArtifactError(
                        "Predicted_ORFs.fasta line {} is not a valid protein sequence".format(
                            line_number
                        )
                    )
                sequence.append(line)
    if header is not None and not sequence:
        raise ArtifactError("Predicted_ORFs.fasta contains an empty record")
    if header is not None:
        fasta_records.append((header, "".join(sequence)))
    if len(fasta_records) != len(expected_records):
        raise ArtifactError(
            "Predicted_ORFs.fasta has {} records; expected {}".format(
                len(fasta_records), len(expected_records)
            )
        )

    expected = {
        (record["locus"], record["strand"]): record
        for record in expected_records
    }
    seen = set()
    for record_number, (raw_header, protein) in enumerate(fasta_records, 1):
        match = FASTA_HEADER.fullmatch(raw_header)
        if not match:
            raise ArtifactError(
                "Predicted_ORFs.fasta record {} has an invalid REPARATION "
                "header".format(record_number)
            )
        identity = (match.group(1), match.group(3))
        if identity not in expected or identity in seen:
            raise ArtifactError(
                "Predicted_ORFs.fasta record {} does not uniquely match the ORF "
                "table".format(record_number)
            )
        seen.add(identity)
        record = expected[identity]
        header_length = int(match.group(4))
        if (
            match.group(2).upper() != record["codon"]
            or header_length != record["length"]
            or len(protein) != record["length"] // 3
        ):
            raise ArtifactError(
                "Predicted_ORFs.fasta record {} has a codon or length "
                "inconsistent with the ORF table".format(record_number)
            )
    if seen != set(expected):
        raise ArtifactError("Predicted_ORFs.fasta records do not match the ORF table")


def _validate_pdf(path):
    _require_file(path, "PDF")
    data = path.read_bytes()
    trailer = re.search(rb"startxref\s+(\d+)\s+%%EOF\s*$", data)
    if not re.match(rb"%PDF-\d\.\d", data) or trailer is None:
        raise ArtifactError("{} is not a valid PDF".format(path.name))
    xref_offset = int(trailer.group(1))
    xref_and_trailer = data[xref_offset : trailer.start()]
    if (
        xref_offset <= 0
        or not xref_and_trailer.startswith(b"xref")
        or b"trailer" not in xref_and_trailer
        or re.search(rb"/Root\s+\d+\s+\d+\s+R", xref_and_trailer) is None
    ):
        raise ArtifactError("{} is not a valid PDF".format(path.name))


def _validate_offsets(path):
    _require_file(path, "P-site offset table")
    with path.open(newline="", encoding="utf-8") as handle:
        numbered_rows = list(enumerate(csv.reader(handle, delimiter="\t"), 1))

    header_index = 0
    while header_index < len(numbered_rows):
        _, row = numbered_rows[header_index]
        if row and not row[0].startswith("#"):
            break
        header_index += 1
    if (
        header_index == len(numbered_rows)
        or numbered_rows[header_index][1] != OFFSET_HEADER
    ):
        raise ArtifactError("p_site_offsets.txt does not have the exact two-column header")
    data_rows = numbered_rows[header_index + 1 :]
    if not data_rows:
        raise ArtifactError("p_site_offsets.txt has no offset data rows")

    lengths = set()
    for line_number, row in data_rows:
        if len(row) != 2:
            raise ArtifactError(
                "p_site_offsets.txt line {} does not have two columns".format(line_number)
            )
        if row[0] == "default":
            length = None
        elif INTEGER.fullmatch(row[0]) and int(row[0]) > 0:
            length = int(row[0])
        else:
            raise ArtifactError(
                "p_site_offsets.txt line {} has an invalid read length".format(line_number)
            )
        if not INTEGER.fullmatch(row[1]) or int(row[1]) < 0:
            raise ArtifactError(
                "p_site_offsets.txt line {} has an invalid P-site offset".format(line_number)
            )
        if length is not None and int(row[1]) >= length:
            raise ArtifactError(
                "p_site_offsets.txt line {} has an offset outside the read length".format(
                    line_number
                )
            )
        if row[0] in lengths:
            raise ArtifactError(
                "p_site_offsets.txt has duplicate length {}".format(row[0])
            )
        lengths.add(row[0])


def _validate_png(path):
    data = path.read_bytes()
    if not data:
        path.unlink()
        print(
            "WARNING: REPARATION produced an empty p_site_offset.png fallback; "
            "it was removed. P-site offsets remain in p_site_offsets.txt.",
            file=sys.stderr,
        )
        return
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ArtifactError("p_site_offset.png is not a valid PNG")

    position = 8
    first_chunk = True
    saw_iend = False
    while position < len(data):
        if position + 12 > len(data):
            raise ArtifactError("p_site_offset.png has a truncated PNG chunk")
        length = int.from_bytes(data[position : position + 4], "big")
        chunk_type = data[position + 4 : position + 8]
        chunk_end = position + 12 + length
        if chunk_end > len(data):
            raise ArtifactError("p_site_offset.png has a truncated PNG chunk")
        chunk_data = data[position + 8 : position + 8 + length]
        stored_crc = int.from_bytes(data[position + 8 + length : chunk_end], "big")
        actual_crc = zlib.crc32(chunk_type + chunk_data) & 0xFFFFFFFF
        if stored_crc != actual_crc:
            raise ArtifactError("p_site_offset.png has an invalid PNG checksum")
        if first_chunk and (chunk_type != b"IHDR" or length != 13):
            raise ArtifactError("p_site_offset.png does not start with a valid IHDR chunk")
        first_chunk = False
        position = chunk_end
        if chunk_type == b"IEND":
            if length != 0 or position != len(data):
                raise ArtifactError("p_site_offset.png has an invalid IEND chunk")
            saw_iend = True
            break
    if not saw_iend:
        raise ArtifactError("p_site_offset.png has no IEND chunk")


def _validate_regular_result_tree(result_dir):
    try:
        root_mode = result_dir.lstat().st_mode
    except FileNotFoundError:
        raise ArtifactError("REPARATION did not produce its result directory")
    if result_dir.is_symlink() or not stat.S_ISDIR(root_mode):
        raise ArtifactError("REPARATION result path is not a real directory")

    for raw_root, directory_names, file_names in os.walk(
        str(result_dir), followlinks=False
    ):
        root = Path(raw_root)
        for name in directory_names:
            path = root / name
            if root == result_dir and name in (
                COMPLETION_RECEIPT,
                COMPLETION_RECEIPT_TEMPORARY,
            ):
                raise ArtifactError(
                    "REPARATION result contains reserved publication path: {}".format(
                        name
                    )
                )
            mode = path.lstat().st_mode
            if path.is_symlink() or not stat.S_ISDIR(mode):
                raise ArtifactError(
                    "REPARATION result contains a symlink or special directory: "
                    "{}".format(path.relative_to(result_dir))
                )
        for name in file_names:
            path = root / name
            if root == result_dir and name in (
                COMPLETION_RECEIPT,
                COMPLETION_RECEIPT_TEMPORARY,
            ):
                raise ArtifactError(
                    "REPARATION result contains reserved publication file: {}".format(
                        name
                    )
                )
            mode = path.lstat().st_mode
            if path.is_symlink() or not stat.S_ISREG(mode):
                raise ArtifactError(
                    "REPARATION result contains a symlink or special file: "
                    "{}".format(path.relative_to(result_dir))
                )


def _validate_result(result_dir, reference_lengths):
    _validate_regular_result_tree(result_dir)
    orf_path = result_dir / "Predicted_ORFs.txt"
    records = _validate_orf_table(orf_path, reference_lengths)
    _validate_bed(result_dir / "Predicted_ORFs.bed", records)
    _validate_fasta(result_dir / "Predicted_ORFs.fasta", records)
    for name in PDF_NAMES:
        _validate_pdf(result_dir / name)
    _validate_offsets(result_dir / "p_site_offsets.txt")
    png = result_dir / "p_site_offset.png"
    if png.exists():
        _validate_png(png)


def _path_exists(path):
    return os.path.lexists(str(path))


def _require_real_directory(path, label):
    if not _path_exists(path):
        return False
    mode = path.lstat().st_mode
    if path.is_symlink() or not stat.S_ISDIR(mode):
        raise ArtifactError("{} is not a real directory: {}".format(label, path))
    return True


def _require_regular_file(path, label):
    if not _path_exists(path):
        return False
    mode = path.lstat().st_mode
    if path.is_symlink() or not stat.S_ISREG(mode):
        raise ArtifactError("{} is not a regular file: {}".format(label, path))
    return True


def _transaction_paths(output_dir):
    prefix = ".{}".format(output_dir.name)
    return {
        "candidate": output_dir.parent / "{}_candidate".format(prefix),
        "backup": output_dir.parent / "{}_backup".format(prefix),
        "marker": output_dir.parent / "{}_transaction".format(prefix),
        "marker_temporary": output_dir.parent
        / "{}_transaction.tmp".format(prefix),
        "lock": output_dir.parent / "{}_lock".format(prefix),
        "staging_locator": output_dir.parent / "{}_staging".format(prefix),
        "staging_locator_temporary": output_dir.parent
        / "{}_staging.tmp".format(prefix),
    }


def _completion_paths(output_dir):
    return {
        "receipt": output_dir / COMPLETION_RECEIPT,
        "temporary": output_dir / COMPLETION_RECEIPT_TEMPORARY,
    }


def _sync_directory(path):
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    descriptor = os.open(str(path), flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _require_private_directory(path, label):
    if not _require_real_directory(path, label):
        return False
    metadata = path.lstat()
    if metadata.st_uid != os.getuid() or stat.S_IMODE(metadata.st_mode) != 0o700:
        raise ArtifactError(
            "{} must be owned by uid {} with mode 0700: {}".format(
                label, os.getuid(), path
            )
        )
    return True


def _require_private_regular_file(path, label):
    if not _require_regular_file(path, label):
        return False
    metadata = path.lstat()
    if metadata.st_uid != os.getuid() or stat.S_IMODE(metadata.st_mode) != 0o600:
        raise ArtifactError(
            "{} must be owned by uid {} with mode 0600: {}".format(
                label, os.getuid(), path
            )
        )
    return True


def _ensure_staging_namespace(temporary_base):
    namespace = _staging_namespace(temporary_base)
    try:
        os.mkdir(str(namespace), 0o700)
    except FileExistsError:
        pass
    _require_private_directory(namespace, "REPARATION staging namespace")
    return namespace


def _write_all(descriptor, content):
    written = 0
    while written < len(content):
        count = os.write(descriptor, content[written:])
        if count == 0:
            raise OSError("short write while creating REPARATION staging metadata")
        written += count


def _write_staging_marker(stage_root, output_dir):
    marker = stage_root / STAGING_MARKER
    content = _staging_marker_content(output_dir)
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(str(marker), flags, 0o600)
    try:
        os.fchmod(descriptor, 0o600)
        _write_all(descriptor, content)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    _sync_directory(stage_root)


def _validate_staging_path(stage_root, output_dir):
    expected_digest = _staging_digest(output_dir)
    raw_path = str(stage_root)
    if (
        not stage_root.is_absolute()
        or not SAFE_PATH.fullmatch(raw_path)
        or ".." in stage_root.parts
        or stage_root.name != "output_{}".format(expected_digest)
        or stage_root.parent.name != "hribo_reparation_{}".format(os.getuid())
    ):
        raise ArtifactError(
            "REPARATION staging locator does not identify this output: {}".format(
                stage_root
            )
        )


def _validate_staging_tree(stage_root, output_dir, locator_authorized=False):
    _validate_staging_path(stage_root, output_dir)
    temporary_base = stage_root.parent.parent
    _require_real_directory(temporary_base, "REPARATION temporary base")
    temporary_base_mode = temporary_base.lstat().st_mode
    if temporary_base_mode & 0o022 and not temporary_base_mode & stat.S_ISVTX:
        raise ArtifactError(
            "REPARATION temporary base is writable and not sticky: {}".format(
                temporary_base
            )
        )
    _require_private_directory(stage_root.parent, "REPARATION staging namespace")
    if not _require_private_directory(stage_root, "REPARATION staging directory"):
        return False

    if locator_authorized:
        # The owned, output-specific locator is durable before mkdir. It is
        # therefore sufficient authority to remove even a partially created or
        # partially deleted tree after an uncatchable process-group SIGKILL.
        return True

    marker = stage_root / STAGING_MARKER
    if not _path_exists(marker):
        # The locator is made durable before mkdir and the marker is the first
        # entry created. An empty tree is therefore the one legitimate
        # pre-marker crash state.
        if any(stage_root.iterdir()):
            raise ArtifactError(
                "REPARATION staging directory has no identity marker: {}".format(
                    stage_root
                )
            )
        return True

    _require_private_regular_file(marker, "REPARATION staging identity marker")
    expected = _staging_marker_content(output_dir)
    metadata = marker.lstat()
    if metadata.st_size != len(expected):
        raise ArtifactError(
            "REPARATION staging identity marker is invalid: {}".format(marker)
        )
    if marker.read_bytes() != expected:
        raise ArtifactError(
            "REPARATION staging identity marker belongs to another output: {}".format(
                marker
            )
        )
    return True


def _remove_staging_tree(stage_root, output_dir, locator_authorized=False):
    if _validate_staging_tree(stage_root, output_dir, locator_authorized):
        shutil.rmtree(str(stage_root))
        _sync_directory(stage_root.parent)


def _staging_locator_content(output_dir, stage_root):
    return "{}\n{}\n{}\n".format(
        STAGING_LOCATOR_VERSION, _staging_digest(output_dir), stage_root
    ).encode("ascii")


def _read_staging_locator(locator, output_dir):
    _require_private_regular_file(locator, "REPARATION staging locator")
    metadata = locator.lstat()
    if metadata.st_size > 8192:
        raise ArtifactError("REPARATION staging locator is invalid: {}".format(locator))
    fields = locator.read_text(encoding="ascii").splitlines()
    if (
        len(fields) != 3
        or fields[0] != STAGING_LOCATOR_VERSION
        or fields[1] != _staging_digest(output_dir)
    ):
        raise ArtifactError("REPARATION staging locator is malformed: {}".format(locator))
    stage_root = Path(fields[2])
    _validate_staging_path(stage_root, output_dir)
    return stage_root


def _write_staging_locator(output_dir, stage_root):
    paths = _transaction_paths(output_dir)
    locator = paths["staging_locator"]
    temporary = paths["staging_locator_temporary"]
    for path, label in (
        (locator, "REPARATION staging locator"),
        (temporary, "temporary REPARATION staging locator"),
    ):
        _require_regular_file(path, label)
        if _path_exists(path):
            raise ArtifactError("{} already exists: {}".format(label, path))

    content = _staging_locator_content(output_dir, stage_root)
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(str(temporary), flags, 0o600)
    try:
        os.fchmod(descriptor, 0o600)
        _write_all(descriptor, content)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(str(temporary), str(locator))
    # This is the durable recovery handle. No staging directory is created
    # until its publication in the output parent has reached disk.
    _sync_directory(output_dir.parent)


def _remove_staging_locator(output_dir, expected_stage_root):
    paths = _transaction_paths(output_dir)
    locator = paths["staging_locator"]
    if _path_exists(locator):
        recorded_stage_root = _read_staging_locator(locator, output_dir)
        if recorded_stage_root != expected_stage_root:
            raise ArtifactError(
                "REPARATION staging locator changed during execution: {}".format(
                    locator
                )
            )
        locator.unlink()
        _sync_directory(output_dir.parent)


def _recover_staging(output_dir):
    paths = _transaction_paths(output_dir)
    temporary = paths["staging_locator_temporary"]
    _require_private_regular_file(
        temporary, "temporary REPARATION staging locator"
    )
    if _path_exists(temporary):
        # Staging data is never created until the final locator is durable.
        temporary.unlink()
        _sync_directory(output_dir.parent)

    locator = paths["staging_locator"]
    _require_regular_file(locator, "REPARATION staging locator")
    if _path_exists(locator):
        old_stage_root = _read_staging_locator(locator, output_dir)
        if _path_exists(old_stage_root):
            _remove_staging_tree(
                old_stage_root, output_dir, locator_authorized=True
            )
        locator.unlink()
        _sync_directory(output_dir.parent)

    temporary_base = _safe_temporary_base()
    _ensure_staging_namespace(temporary_base)
    stage_root = _staging_path(output_dir, temporary_base)
    # This fallback also recovers a directory whose locator was manually lost,
    # provided its output-specific marker is intact (or mkdir was interrupted
    # before the marker could be written).
    if _path_exists(stage_root):
        _remove_staging_tree(stage_root, output_dir)
    return stage_root


@contextmanager
def _staging_directory(output_dir):
    stage_root = _recover_staging(output_dir)
    _write_staging_locator(output_dir, stage_root)
    try:
        os.mkdir(str(stage_root), 0o700)
        _require_private_directory(stage_root, "REPARATION staging directory")
        _write_staging_marker(stage_root, output_dir)
        yield stage_root
    finally:
        # SIGKILL skips this block, but the durable locator above lets the next
        # lock owner perform the identical bounded cleanup.
        if _path_exists(stage_root):
            _remove_staging_tree(stage_root, output_dir, locator_authorized=True)
        _remove_staging_locator(output_dir, stage_root)


def _invalidate_completion_receipt(output_dir):
    if not _require_real_directory(output_dir, "published output"):
        return
    paths = _completion_paths(output_dir)
    for path, label in (
        (paths["temporary"], "temporary completion receipt"),
        (paths["receipt"], "completion receipt"),
    ):
        _require_regular_file(path, label)
    removed = False
    for path in (paths["temporary"], paths["receipt"]):
        if _path_exists(path):
            path.unlink()
            removed = True
    if removed:
        # This is the durable barrier before any result-tree mutation: after it
        # returns, a crash cannot leave a completion signal for an in-progress
        # or failed rerun.
        _sync_directory(output_dir)


def _publish_completion_receipt(output_dir):
    if not _require_real_directory(output_dir, "published output"):
        raise ArtifactError(
            "cannot publish a completion receipt without a result directory"
        )
    paths = _completion_paths(output_dir)
    if _path_exists(paths["receipt"]) or _path_exists(paths["temporary"]):
        raise ArtifactError("REPARATION completion receipt already exists")

    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(str(paths["temporary"]), flags, 0o600)
    try:
        written = 0
        while written < len(COMPLETION_RECEIPT_CONTENT):
            count = os.write(
                descriptor, COMPLETION_RECEIPT_CONTENT[written:]
            )
            if count == 0:
                raise OSError("short write while creating completion receipt")
            written += count
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(str(paths["temporary"]), str(paths["receipt"]))
    # The result tree is already committed and synced. Persist the final rename
    # separately so the receipt can never precede the artifacts it represents.
    _sync_directory(output_dir)


def _sync_tree(path):
    for raw_root, _, file_names in os.walk(str(path), topdown=False):
        root = Path(raw_root)
        for name in file_names:
            file_path = root / name
            flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
            descriptor = os.open(str(file_path), flags)
            try:
                if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                    raise ArtifactError(
                        "publication candidate contains a special file: {}".format(
                            file_path.relative_to(path)
                        )
                    )
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        _sync_directory(root)


def _remove_transaction_directory(path, label):
    if _require_real_directory(path, label):
        shutil.rmtree(str(path))


def _legacy_transaction_paths(output_dir):
    backup_prefix = ".{}_backup_".format(output_dir.name)
    publish_prefix = ".{}_publish_".format(output_dir.name)
    backups = []
    publish_roots = []
    for path in output_dir.parent.iterdir():
        if path.name.startswith(backup_prefix):
            _require_real_directory(path, "legacy publication backup")
            backups.append(path)
        elif path.name.startswith(publish_prefix):
            _require_real_directory(path, "legacy publication candidate")
            publish_roots.append(path)
    return backups, publish_roots


def _recover_legacy_publication(output_dir):
    backups, publish_roots = _legacy_transaction_paths(output_dir)
    if not _path_exists(output_dir) and len(backups) > 1:
        raise ArtifactError(
            "multiple legacy publication backups exist for {}; refusing to guess".format(
                output_dir
            )
        )
    if not _path_exists(output_dir) and len(backups) == 1:
        os.replace(str(backups[0]), str(output_dir))
        backups = []
        _sync_directory(output_dir.parent)
    for path in backups:
        _remove_transaction_directory(path, "legacy publication backup")
    for path in publish_roots:
        _remove_transaction_directory(path, "legacy publication candidate")
    if backups or publish_roots:
        _sync_directory(output_dir.parent)


def _read_marker(path):
    _require_regular_file(path, "publication transaction marker")
    state = path.read_text(encoding="ascii")
    if state not in ("present\n", "absent\n"):
        raise ArtifactError("publication transaction marker is malformed: {}".format(path))
    return state.rstrip("\n")


def _recover_publication(output_dir):
    paths = _transaction_paths(output_dir)
    candidate = paths["candidate"]
    backup = paths["backup"]
    marker = paths["marker"]
    marker_temporary = paths["marker_temporary"]

    _require_real_directory(candidate, "publication candidate")
    _require_real_directory(backup, "publication backup")
    _require_regular_file(marker_temporary, "temporary publication marker")

    if _path_exists(marker):
        state = _read_marker(marker)
        if state == "present":
            if _path_exists(backup):
                if _path_exists(output_dir):
                    _require_real_directory(output_dir, "partially published output")
                    shutil.rmtree(str(output_dir))
                os.replace(str(backup), str(output_dir))
                _sync_directory(output_dir.parent)
            elif not _path_exists(output_dir):
                raise ArtifactError(
                    "interrupted publication has neither its original output nor backup"
                )
        else:
            if _path_exists(backup):
                raise ArtifactError(
                    "first publication unexpectedly contains an original-output backup"
                )
            if _path_exists(output_dir):
                _require_real_directory(output_dir, "partially published output")
                shutil.rmtree(str(output_dir))
                _sync_directory(output_dir.parent)
        _remove_transaction_directory(candidate, "publication candidate")
        marker.unlink()
        _sync_directory(output_dir.parent)
    elif _path_exists(backup):
        if not _path_exists(output_dir):
            raise ArtifactError(
                "committed publication backup exists but the published output is missing"
            )
        _require_real_directory(output_dir, "published output")
        _remove_transaction_directory(backup, "publication backup")
        _remove_transaction_directory(candidate, "publication candidate")
        _sync_directory(output_dir.parent)
    else:
        _remove_transaction_directory(candidate, "publication candidate")

    if _path_exists(marker_temporary):
        marker_temporary.unlink()
        _sync_directory(output_dir.parent)
    _recover_legacy_publication(output_dir)


@contextmanager
def _publication_lock(output_dir):
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    if output_dir.parent.is_symlink():
        raise ArtifactError(
            "refusing to publish through a symlinked output parent: {}".format(
                output_dir.parent
            )
        )
    lock_path = _transaction_paths(output_dir)["lock"]
    flags = os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(str(lock_path), flags, 0o600)
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ArtifactError(
                "publication lock is not a regular file: {}".format(lock_path)
            )
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        yield
    finally:
        os.close(descriptor)


def _write_marker(paths, state):
    marker_temporary = paths["marker_temporary"]
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(str(marker_temporary), flags, 0o600)
    try:
        os.write(descriptor, (state + "\n").encode("ascii"))
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(str(marker_temporary), str(paths["marker"]))
    _sync_directory(paths["marker"].parent)


def _publish(staged_result, output_dir, reference_lengths):
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    paths = _transaction_paths(output_dir)
    candidate = paths["candidate"]
    backup = paths["backup"]
    marker_written = False
    try:
        shutil.copytree(str(staged_result), str(candidate))
        _validate_result(candidate, reference_lengths)
        _sync_tree(candidate)
        _sync_directory(output_dir.parent)
        old_output = _path_exists(output_dir)
        _write_marker(paths, "present" if old_output else "absent")
        marker_written = True
        if old_output:
            os.replace(str(output_dir), str(backup))
            _sync_directory(output_dir.parent)
        os.replace(str(candidate), str(output_dir))
        _sync_directory(output_dir.parent)
        paths["marker"].unlink()
        marker_written = False
    except BaseException:
        if marker_written or _path_exists(paths["marker"]):
            _recover_publication(output_dir)
        else:
            _remove_transaction_directory(candidate, "publication candidate")
        raise

    # Removing the marker is the result-tree commit point. Leave any backup for
    # deterministic startup reconciliation and keep the complete new tree; the
    # caller publishes the separate completion receipt only after this returns.
    try:
        _sync_directory(output_dir.parent)
        _remove_transaction_directory(backup, "publication backup")
        _sync_directory(output_dir.parent)
    except Exception as error:
        print(
            "WARNING: REPARATION output was committed, but transaction cleanup "
            "will be retried on the next run: {}".format(error),
            file=sys.stderr,
        )


def _run_engine(args, stage_root, aliases):
    result_dir = stage_root / "result"
    engine = args.engine
    if os.sep in engine:
        engine = str(Path(engine).expanduser().resolve())
    command = [
        engine,
        "-bam",
        str(aliases["bam"]),
        "-g",
        str(aliases["genome"]),
        "-gtf",
        str(aliases["gtf"]),
        "-db",
        str(aliases["database"]),
        "-wdir",
        str(result_dir),
        "-threads",
        str(args.threads),
    ]
    completed = subprocess.run(command, cwd=str(stage_root))
    return completed.returncode, result_dir


def main():
    args = _arguments()
    if args.threads < 1:
        print("ERROR: --threads must be at least 1", file=sys.stderr)
        return 2

    try:
        genome = _input_file(args.genome, "genome")
        gtf = _input_file(args.gtf, "GTF")
        database = _input_file(args.database, "protein database")
        bam = _input_file(args.bam, "BAM")
        bai = _input_file(args.bai, "BAM index")
        output_dir = _output_directory(args.output_dir)
        _validate_output_inputs(
            output_dir,
            {
                "genome": genome,
                "GTF": gtf,
                "protein database": database,
                "BAM": bam,
                "BAM index": bai,
            },
        )
        with _publication_lock(output_dir):
            # Invalidate both before and after recovery. A pre-commit backup
            # created by an older wrapper may itself contain a receipt, so the
            # second pass guarantees the engine never starts behind a stale
            # completion signal.
            _invalidate_completion_receipt(output_dir)
            _recover_publication(output_dir)
            _invalidate_completion_receipt(output_dir)
            reference_lengths = _reference_lengths(genome)
            with _staging_directory(output_dir) as stage_root:
                if not SAFE_PATH.fullmatch(str(stage_root)):
                    raise ArtifactError(
                        "temporary staging path is not safe for the legacy REPARATION engine"
                    )
                aliases = _stage_inputs(stage_root, genome, gtf, database, bam, bai)
                returncode, result_dir = _run_engine(args, stage_root, aliases)
                if returncode != 0:
                    return returncode
                _validate_result(result_dir, reference_lengths)
                _publish(result_dir, output_dir, reference_lengths)
            # Staging cleanup completes before the separate success signal, so
            # a cleanup failure can never leave a misleading receipt.
            _publish_completion_receipt(output_dir)
    except (ArtifactError, OSError, shutil.Error, UnicodeError, csv.Error) as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
