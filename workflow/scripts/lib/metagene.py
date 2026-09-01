"""
Contains scripts for metagene profiling analysis.
Author: Rick Gelhausen
"""

import lib.misc as misc
import numpy as np


def _transcript_index(genomic_position, window_start, window_stop, strand):
    """Return an array index that always increases in the 5'-to-3' direction."""
    if strand == "+":
        return genomic_position - window_start
    return window_stop - genomic_position


def _window_bounds(codon, strand, anchor, positions_out_ORF, positions_in_ORF):
    """Build an inclusive genomic window around a start or stop codon."""
    start, stop = codon

    if anchor == "start":
        if strand == "+":
            return start - positions_out_ORF, start + positions_in_ORF - 1
        return stop - positions_in_ORF + 1, stop + positions_out_ORF

    if strand == "+":
        return stop - positions_in_ORF + 1, stop + positions_out_ORF
    return start - positions_out_ORF, start + positions_in_ORF - 1


def _point_position(read_interval, strand, mapping_method):
    """Select the requested biological read position in genomic coordinates."""
    if mapping_method == "fiveprime":
        return read_interval[0] if strand == "+" else read_interval[1]
    if mapping_method == "threeprime":
        return read_interval[1] if strand == "+" else read_interval[0]
    if mapping_method == "centered":
        return round((read_interval[0] + read_interval[1]) / 2)
    return None


def _metagene_mapping(
    codon_dict,
    read_interval_dict,
    positions_out_ORF,
    positions_in_ORF,
    mapping_method,
    anchor,
):
    """Map reads around one codon type onto a transcript-oriented profile."""
    coverage_mapping = {}
    window_length = positions_out_ORF + positions_in_ORF

    for strand in codon_dict:
        for chromosome in codon_dict[strand]:
            chromosome_coverage = coverage_mapping.setdefault(chromosome, {})

            for codon in codon_dict[strand][chromosome]:
                window_start, window_stop = _window_bounds(
                    codon,
                    strand,
                    anchor,
                    positions_out_ORF,
                    positions_in_ORF,
                )

                read_index = read_interval_dict.get((chromosome, strand))
                if read_index is None:
                    # An annotation may retain a contig/strand on which no read
                    # survived filtering.  That is valid zero evidence, not a
                    # malformed interval dictionary.
                    continue

                reads = read_index.find((window_start, window_stop))
                for read_interval in reads:
                    read_length = read_interval[2]
                    if read_length not in chromosome_coverage:
                        chromosome_coverage[read_length] = np.zeros(
                            window_length, dtype=np.intp
                        )
                    coverage = chromosome_coverage[read_length]

                    if mapping_method == "global":
                        intersection = misc.get_overlap_bounderies(
                            (window_start, window_stop), read_interval
                        )
                        first_index = _transcript_index(
                            intersection[0], window_start, window_stop, strand
                        )
                        last_index = _transcript_index(
                            intersection[1], window_start, window_stop, strand
                        )
                        low_index, high_index = sorted((first_index, last_index))
                        coverage[low_index : high_index + 1] += 1
                        continue

                    point = _point_position(read_interval, strand, mapping_method)
                    if point is not None and window_start <= point <= window_stop:
                        index = _transcript_index(
                            point, window_start, window_stop, strand
                        )
                        coverage[index] += 1

    return coverage_mapping


def metagene_mapping_start(
    start_codon_dict,
    read_interval_dict,
    positions_out_ORF,
    positions_in_ORF,
    mapping_method,
):
    """Map reads around starts onto ``[-outside, inside)`` coordinates."""
    return _metagene_mapping(
        start_codon_dict,
        read_interval_dict,
        positions_out_ORF,
        positions_in_ORF,
        mapping_method,
        "start",
    )


def metagene_mapping_stop(
    stop_codon_dict,
    read_interval_dict,
    positions_out_ORF,
    positions_in_ORF,
    mapping_method,
):
    """Map reads around stops onto ``[-inside, outside)`` coordinates."""
    return _metagene_mapping(
        stop_codon_dict,
        read_interval_dict,
        positions_out_ORF,
        positions_in_ORF,
        mapping_method,
        "stop",
    )
