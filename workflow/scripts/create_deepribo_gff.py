#!/usr/bin/env python3
"""Convert DeepRibo predictions to reference-validated GFF3 records."""

import argparse
import collections
import csv
import io
import math
import os
import re
import tempfile
from decimal import Decimal, InvalidOperation
from pathlib import Path

import pandas as pd

import gff_utils


GFF3_HEADER = "##gff-version 3\n"
GFF_COLUMNS = [
    "seqName",
    "source",
    "type",
    "start",
    "stop",
    "score",
    "strand",
    "phase",
    "attribute",
]
INTERVAL = re.compile(r"^(\d+)-(\d+)$")
GENOME_SEQUENCE = re.compile(r"^[ACGTRYSWKMBDHVNU]+$", re.IGNORECASE)
REQUIRED_COLUMNS = {"strand", "locus", "pred", "dist", "SS_pred_rank"}


class PredictionError(ValueError):
    """A DeepRibo prediction cannot be represented against the reference."""


def fasta_lengths(path):
    """Return FASTA lengths under DeepRibo's first-token ID convention."""

    lengths = {}
    current = None
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields:
                    raise PredictionError(
                        f"{path}: line {line_number} has no FASTA identifier"
                    )
                current = fields[0]
                if current in lengths:
                    raise PredictionError(
                        f"{path}: duplicate FASTA identifier {current!r}"
                    )
                lengths[current] = 0
            elif current is None:
                raise PredictionError(
                    f"{path}: sequence data precedes the first FASTA header"
                )
            else:
                if not GENOME_SEQUENCE.fullmatch(line):
                    raise PredictionError(
                        f"{path}: line {line_number} contains non-IUPAC "
                        "FASTA sequence data"
                    )
                lengths[current] += len(line)

    if not lengths:
        raise PredictionError(f"{path}: no FASTA records found")
    empty = [identifier for identifier, length in lengths.items() if length == 0]
    if empty:
        raise PredictionError(f"{path}: FASTA record {empty[0]!r} has no sequence")
    return lengths


def parse_locus(locus, line_number):
    """Parse ``contig:start-stop`` while permitting colons in the contig ID."""

    if ":" not in locus:
        raise PredictionError(
            f"predictions.csv line {line_number} has invalid locus {locus!r}"
        )
    chromosome, interval = locus.rsplit(":", 1)
    match = INTERVAL.fullmatch(interval)
    if not chromosome or match is None:
        raise PredictionError(
            f"predictions.csv line {line_number} has invalid locus {locus!r}"
        )
    start, stop = int(match.group(1)), int(match.group(2))
    if start < 1 or stop < start:
        raise PredictionError(
            f"predictions.csv line {line_number} has invalid coordinate range "
            f"{start}-{stop}"
        )
    return chromosome, start, stop


def complete_orf_coordinates(chromosome, start, stop, strand, lengths, line_number):
    """Restore DeepRibo's two omitted stop-codon bases and check the result."""

    if strand == "+":
        complete_start, complete_stop = start, stop + 2
    elif strand == "-":
        complete_start, complete_stop = start - 2, stop
    else:
        raise PredictionError(
            f"predictions.csv line {line_number} has invalid strand {strand!r}"
        )

    genome_length = lengths.get(chromosome)
    if genome_length is None:
        raise PredictionError(
            f"predictions.csv line {line_number} names unknown reference "
            f"sequence {chromosome!r}"
        )
    if complete_start < 1 or complete_stop > genome_length:
        raise PredictionError(
            f"predictions.csv line {line_number} complete ORF "
            f"{chromosome}:{complete_start}-{complete_stop}:{strand} is outside "
            f"the reference sequence length {genome_length}"
        )

    complete_length = complete_stop - complete_start + 1
    if complete_length < 6:
        raise PredictionError(
            f"predictions.csv line {line_number} complete ORF "
            f"{chromosome}:{complete_start}-{complete_stop}:{strand} has length "
            f"{complete_length}; a complete CDS must be at least 6 nt"
        )
    if complete_length % 3 != 0:
        raise PredictionError(
            f"predictions.csv line {line_number} complete ORF "
            f"{chromosome}:{complete_start}-{complete_stop}:{strand} has length "
            f"{complete_length}, which is not divisible by 3"
        )
    return complete_start, complete_stop


def finite_number(value, field, line_number):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise PredictionError(
            f"predictions.csv line {line_number} has non-numeric {field}"
        )
    if not math.isfinite(number):
        raise PredictionError(
            f"predictions.csv line {line_number} has non-finite {field}"
        )
    return number


def integer_number(value, field, line_number):
    try:
        number = Decimal(str(value).strip())
    except (InvalidOperation, ValueError):
        raise PredictionError(
            f"predictions.csv line {line_number} has non-numeric {field}"
        )
    if not number.is_finite():
        raise PredictionError(
            f"predictions.csv line {line_number} has non-finite {field}"
        )
    if number != number.to_integral_value():
        raise PredictionError(
            f"predictions.csv line {line_number} has non-integer {field}"
        )
    return int(number)


def to_gff3(args):
    # Keep numeric metadata as text until it has been validated.  Letting pandas
    # infer float64 here silently rounds integer-valued decimals above 2**53.
    input_df = pd.read_csv(
        args.predictedORFs, sep=",", dtype=str, keep_default_na=False
    )
    missing = sorted(REQUIRED_COLUMNS - set(input_df.columns))
    if missing:
        raise PredictionError(
            "predictions.csv is missing required column(s): " + ", ".join(missing)
        )
    reference_lengths = fasta_lengths(args.genome)
    record = collections.namedtuple("Pandas", GFF_COLUMNS)

    rows = []
    for line_number, row in enumerate(
        input_df.itertuples(index=False, name="Pandas"), start=2
    ):
        strand = str(getattr(row, "strand"))
        locus = str(getattr(row, "locus"))
        pred = finite_number(getattr(row, "pred"), "pred", line_number)
        dist = integer_number(getattr(row, "dist"), "dist", line_number)
        ss_pred_rank = integer_number(
            getattr(row, "SS_pred_rank"), "SS_pred_rank", line_number
        )
        if ss_pred_rank < 0:
            raise PredictionError(
                f"predictions.csv line {line_number} has negative SS_pred_rank"
            )

        chromosome, start, stop = parse_locus(locus, line_number)
        start, stop = complete_orf_coordinates(
            chromosome,
            start,
            stop,
            strand,
            reference_lengths,
            line_number,
        )

        if ss_pred_rank == 999999:
            continue

        identifier = f"{chromosome}:{start}-{stop}:{strand}"
        attribute = gff_utils.format_attributes(
            [
                ("ID", identifier),
                ("pred_value", str(pred)),
                ("deepribo_distance", str(dist)),
                ("method", "deepribo"),
                ("condition", args.condition),
                ("replicate", args.replicate),
            ]
        )

        # DeepRibo's distance is prediction metadata, not the CDS reading frame.
        # Predictions are complete ORFs beginning at their start codon, so their
        # required CDS phase is zero on either strand.
        rows.append(
            record(
                chromosome,
                "deepribo",
                "CDS",
                start,
                stop,
                pred,
                strand,
                "0",
                attribute,
            )
        )

    return pd.DataFrame.from_records(rows, columns=GFF_COLUMNS)


def render_gff(dataframe):
    output = io.StringIO()
    output.write(GFF3_HEADER)
    dataframe.to_csv(
        output,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
    )
    return output.getvalue()


def atomic_write(path, content):
    """Replace the output only after its complete new contents are staged."""

    output = Path(path)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp", text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fchmod(handle.fileno(), 0o644)
            os.fsync(handle.fileno())
        os.replace(temporary_name, output)
        directory_descriptor = os.open(
            str(output.parent), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
        )
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Convert DeepRibo predictions to reference-validated GFF3."
    )
    parser.add_argument(
        "-i",
        "--inputCSV",
        dest="predictedORFs",
        required=True,
        help="DeepRibo predictions.csv",
    )
    parser.add_argument(
        "-g",
        "--genome",
        required=True,
        help="reference genome FASTA used to generate the predictions",
    )
    parser.add_argument(
        "-c", "--condition", required=True, help="condition of the current file"
    )
    parser.add_argument(
        "-r", "--replicate", required=True, help="replicate of the current file"
    )
    parser.add_argument(
        "-o",
        "--outputGFF",
        dest="outputGFF",
        required=True,
        help="output GFF3 file",
    )
    args = parser.parse_args()

    try:
        gff3df = to_gff3(args)
        # Explicit tie-breakers and a stable sort: sorting on score alone left rows
        # of equal score in an order that varied with the pandas build.
        gff3df = gff3df.sort_values(
            by=["score", "seqName", "start", "stop", "strand"], kind="stable"
        )
        atomic_write(args.outputGFF, render_gff(gff3df))
    except (OSError, TypeError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"create_deepribo_gff.py: error: {error}\n")


if __name__ == "__main__":
    main()
