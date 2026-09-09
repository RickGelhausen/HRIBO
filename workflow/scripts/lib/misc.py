"""
Contains miscellaneous scripts for different calculations.
Author: Rick Gelhausen
"""
import csv

import pandas as pd
import numpy as np


def read_mapped_read_summary(path):
    """Parse and aggregate a three-column mapped-read summary.

    ``total_mapped_reads.py`` intentionally writes one row per library and
    contig.  Normalisation denominators, however, describe the complete
    library, so consumers must sum those rows rather than selecting the row for
    the feature's contig.

    Returns ``(counts_by_library_and_contig, totals_by_library)``.  Individual
    contigs may contain zero reads, but every library represented in the file
    must have a positive aggregate count.
    """
    counts_by_contig = {}
    totals_by_library = {}

    with open(path, newline="") as handle:
        for line_number, row in enumerate(
            csv.reader(handle, delimiter="\t"), start=1
        ):
            if not row or all(not field.strip() for field in row):
                continue
            if len(row) != 3:
                raise ValueError(
                    f"{path}: line {line_number}: expected three tab-separated "
                    "fields (library, contig, mapped reads)"
                )

            library, contig, raw_count = (field.strip() for field in row)
            if not library or not contig:
                raise ValueError(
                    f"{path}: line {line_number}: library and contig must not be empty"
                )
            try:
                count = float(raw_count)
            except ValueError as error:
                raise ValueError(
                    f"{path}: line {line_number}: mapped-read count "
                    f"{raw_count!r} is not numeric"
                ) from error
            if not np.isfinite(count):
                raise ValueError(
                    f"{path}: line {line_number}: mapped-read count must be finite"
                )
            if count < 0:
                raise ValueError(
                    f"{path}: line {line_number}: mapped-read count must not be negative"
                )

            key = (library, contig)
            if key in counts_by_contig:
                raise ValueError(
                    f"{path}: line {line_number}: duplicate mapped-read row for "
                    f"library {library!r} and contig {contig!r}"
                )
            counts_by_contig[key] = count
            totals_by_library[library] = totals_by_library.get(library, 0) + count

    if not counts_by_contig:
        raise ValueError(f"{path}: mapped-read summary contains no data rows")

    for library, total in totals_by_library.items():
        if total == 0:
            raise ValueError(
                f"{path}: library {library!r} has zero mapped reads across all contigs; "
                "normalisation is undefined"
            )

    return counts_by_contig, totals_by_library


def library_read_total(counts_by_contig):
    """Return a validated, positive library-wide total from contig counts."""
    if not counts_by_contig:
        raise ValueError(
            "mapped-read counts are missing; library-wide normalisation is undefined"
        )

    total = 0
    for contig, count in counts_by_contig.items():
        if isinstance(count, bool):
            raise ValueError(
                f"mapped-read count for contig {contig!r} must be numeric"
            )
        try:
            numeric_count = float(count)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"mapped-read count for contig {contig!r} must be numeric"
            ) from error
        if not np.isfinite(numeric_count):
            raise ValueError(
                f"mapped-read count for contig {contig!r} must be finite"
            )
        if numeric_count < 0:
            raise ValueError(
                f"mapped-read count for contig {contig!r} must not be negative"
            )
        total += numeric_count

    if total == 0:
        raise ValueError(
            "library has zero mapped reads across all contigs; normalisation is undefined"
        )
    return total

def get_overlap_bounderies(a, b):
    """
    get the overlap between two intervals
    """
    return (max(a[0], b[0]), min(a[1], b[1]))


def get_aligned_blocks(read_interval):
    """Return inclusive CIGAR-aligned blocks from an interval tuple.

    Three-field tuples predate block-aware alignment parsing and remain useful
    in tests and callers that construct synthetic intervals directly.
    """
    if len(read_interval) >= 4:
        return read_interval[3]
    return ((read_interval[0], read_interval[1]),)

def calculate_rpkm(gene_length, read_counts, total_counts):
    """
    Calculate rpkm for a gene
    """

    if gene_length <= 0:
        raise ValueError("gene length must be positive to calculate RPKM")
    if total_counts <= 0:
        raise ValueError(
            "library-wide mapped-read total must be positive to calculate RPKM"
        )
    return (read_counts * 1000000000) / (total_counts * gene_length)

def count_reads(read_intervals_dict, chromosome, strand, beginning, end, mapping_method):
    """
    Select the right reads based on the mapping method.
    """
    read_intervals = read_intervals_dict[(chromosome,strand)].find((beginning, end))
    allowed_reads = []
    for read_interval in read_intervals:
        if mapping_method == "global":
            if any(
                block_start <= end and block_stop >= beginning
                for block_start, block_stop in get_aligned_blocks(read_interval)
            ):
                allowed_reads.append(read_interval)
        elif mapping_method == "threeprime":
            if strand == "+":
                if read_interval[1] <= end:
                    allowed_reads.append(read_interval)
            else:
                if read_interval[0] >= beginning:
                    allowed_reads.append(read_interval)
        elif mapping_method == "fiveprime":
            if strand == "+":
                if read_interval[0] >= beginning:
                    allowed_reads.append(read_interval)
            else:
                if read_interval[1] <= end:
                    allowed_reads.append(read_interval)
        elif mapping_method == "centered":
            center = round((read_interval[0] + read_interval[1]) / 2)
            if center >= beginning and center <= end:
                allowed_reads.append(read_interval)

    return len(allowed_reads)

def normalize_coverage(metagene_coverage_dict, total_counts_dict):
    """
    Per million normalization for the metagene coverage
    """

    library_total = library_read_total(total_counts_dict)
    for chrom in metagene_coverage_dict:
        for read_length in metagene_coverage_dict[chrom]:
            metagene_coverage_dict[chrom][read_length] = (
                metagene_coverage_dict[chrom][read_length] / library_total
            ) * 1000000

    return metagene_coverage_dict


def retain_read_lengths(metagene_coverage_dict, read_lengths):
    """Keep only configured read lengths and discard now-empty contigs."""
    wanted = {int(read_length) for read_length in read_lengths}
    return {
        chromosome: {
            read_length: values
            for read_length, values in coverage.items()
            if int(read_length) in wanted
        }
        for chromosome, coverage in metagene_coverage_dict.items()
        if any(int(read_length) in wanted for read_length in coverage)
    }

def window_normalize_df(df, window_size):
    """
    Normalize every read length column by the total read length number and the window size.
    """

    columns = df.columns[1:].tolist()
    for column in columns:
        scale = df[column].sum() / window_size
        if scale == 0:
            # A read length filled in by ``equalize_dictionary_keys`` represents
            # explicit zero evidence.  Dividing that column by zero would turn a
            # valid empty profile into NaNs in both plots and workbooks.
            df[column] = 0.0
        else:
            df[column] = df[column].div(scale)

    return df

def create_data_frame(
    metagene_dict,
    positions_out_ORF,
    positions_in_ORF,
    state,
    read_lengths=None,
):
    """
    Create a data frame containing the metagene profiling read counts.
    """

    dataframe_dict = {}
    if state == "start":
        coordinates = list(range(-positions_out_ORF, positions_in_ORF, 1))
    else:
        coordinates = list(range(-positions_in_ORF, positions_out_ORF, 1))

    if not metagene_dict:
        # Keep empty scientific evidence explicit in the workbook instead of
        # relying on the spreadsheet engine's anonymous fallback sheet.
        frame = pd.DataFrame({"coordinates": coordinates})
        for read_length in sorted({int(value) for value in read_lengths or []}):
            frame[str(read_length)] = 0
        dataframe_dict["no_evidence"] = frame
        return dataframe_dict

    for chrom in metagene_dict:
        if chrom not in dataframe_dict:
            dataframe_dict[chrom] = pd.DataFrame()
            dataframe_dict[chrom]["coordinates"] = coordinates
        for read_length in sorted(metagene_dict[chrom].keys()):
            dataframe_dict[chrom][f"{read_length}"] = metagene_dict[chrom][read_length]

    return dataframe_dict

def equalize_dictionary_keys(
    start_dict,
    stop_dict,
    positions_out_ORF,
    positions_in_ORF,
    read_lengths=None,
):
    """
    Ensure that both dictionaries have the same set of keys.
    Create new keys for missing values and initialize them with list of 0s.
    """
    window_length = positions_out_ORF + positions_in_ORF

    def empty_window():
        # A fresh array per read length: sharing one array would alias every
        # filled-in read length to the same buffer.
        return np.zeros(window_length, dtype=np.intp)

    unique_keys = set(start_dict.keys()).union(set(stop_dict.keys()))

    for key in unique_keys:
        for coverage_dict in (start_dict, stop_dict):
            if key not in coverage_dict:
                coverage_dict[key] = {}

    if read_lengths is None:
        observed_lengths = {
            int(read_length)
            for coverage_dict in (start_dict, stop_dict)
            for chromosome_coverage in coverage_dict.values()
            for read_length in chromosome_coverage
        }
        if not observed_lengths:
            return start_dict, stop_dict
        balanced_lengths = range(min(observed_lengths), max(observed_lengths) + 1)
    else:
        balanced_lengths = sorted({int(read_length) for read_length in read_lengths})

    for key in unique_keys:
        for coverage_dict in (start_dict, stop_dict):
            for read_length in balanced_lengths:
                if read_length not in coverage_dict[key]:
                    coverage_dict[key][read_length] = empty_window()

    return start_dict, stop_dict
