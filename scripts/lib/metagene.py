"""
Contains scripts for metagene profiling analysis.
Author: Rick Gelhausen
"""

import lib.misc as misc
import numpy as np

def metagene_mapping_start(start_codon_dict, read_interval_dict, positions_out_ORF, positions_in_ORF, mapping_method):
    """
    For every start codon map the read to the positions in the desired interval.
    """

    coverage_mapping = {}
    window_length = positions_out_ORF + positions_in_ORF
    for strand in start_codon_dict:
        for chromosome in start_codon_dict[strand]:
            if chromosome not in coverage_mapping:
                coverage_mapping[chromosome] = {}
            for start_codon in start_codon_dict[strand][chromosome]:
                start, stop = start_codon

                if strand == "+":
                    window_start = start - positions_out_ORF
                    window_stop = start + positions_in_ORF - 1
                    for read_interval in read_interval_dict[(chromosome, strand)].find((window_start, window_stop)):
                        intersection = misc.get_overlap_bounderies((window_start, window_stop), read_interval)

                        read_length = read_interval[2]
                        if read_length not in coverage_mapping[chromosome]:
                            coverage_mapping[chromosome][read_length] = np.zeros(window_length, dtype=np.int)

                        coverage_low_bound = intersection[0] - window_start
                        coverage_high_bound = intersection[1] - window_start

                        if mapping_method == "fiveprime":
                            if read_interval[0] >= window_start:
                                coverage_mapping[chromosome][read_length][coverage_low_bound] += 1
                        elif mapping_method == "threeprime":
                            if read_interval[1] <= window_stop:
                                coverage_mapping[chromosome][read_length][coverage_high_bound] += 1

                        elif mapping_method == "centered":
                            center = round((read_interval[0] + read_interval[1]) / 2)
                            center_coverage = center - window_start
                            if center >= window_start and center <= window_stop:
                                coverage_mapping[chromosome][read_length][center_coverage] += 1

                        elif mapping_method == "global":
                            coverage_mapping[chromosome][read_length][coverage_low_bound:coverage_high_bound+1] += 1

                else:
                    window_start = stop - positions_in_ORF + 1
                    window_stop = stop + positions_out_ORF
                    for read_interval in read_interval_dict[(chromosome, strand)].find((window_start, window_stop)):
                        intersection = misc.get_overlap_bounderies((window_start, window_stop), read_interval)

                        coverage_low_bound = (window_length - 1) - abs(intersection[0] - window_start)
                        coverage_high_bound = (window_length - 1) - abs(intersection[1] - window_start)

                        read_length = read_interval[2]
                        if read_length not in coverage_mapping[chromosome]:
                            coverage_mapping[chromosome][read_length] = np.zeros(window_length, dtype=np.int)

                        if mapping_method == "fiveprime":
                            if read_interval[1] <= window_stop:
                                coverage_mapping[chromosome][read_length][coverage_high_bound] += 1
                        elif mapping_method == "threeprime":
                            if read_interval[0] >= window_start:
                                coverage_mapping[chromosome][read_length][coverage_low_bound] += 1

                        elif mapping_method == "centered":
                            center = round((read_interval[0] + read_interval[1]) / 2)
                            center_coverage = (window_length - 1) - abs(center - window_start)
                            if center >= window_start and center <= window_stop:
                                coverage_mapping[chromosome][read_length][center_coverage] += 1

                        elif mapping_method == "global":
                            coverage_mapping[chromosome][read_length][coverage_low_bound:coverage_high_bound+1] += 1

    return coverage_mapping

def metagene_mapping_stop(stop_codon_dict, read_interval_dict, positions_out_ORF, positions_in_ORF, mapping_method):
    """
    For every start codon map the read to the positions in the desired interval.
    """

    coverage_mapping = {}
    window_length = positions_out_ORF + positions_in_ORF
    for strand in stop_codon_dict:
        for chromosome in stop_codon_dict[strand]:
            if chromosome not in coverage_mapping:
                coverage_mapping[chromosome] = {}
            for stop_codon in stop_codon_dict[strand][chromosome]:
                start, stop = stop_codon

                if strand == "+":
                    window_start = stop - positions_in_ORF + 1
                    window_stop = stop + positions_out_ORF
                    for read_interval in read_interval_dict[(chromosome, strand)].find((window_start, window_stop)):
                        intersection = misc.get_overlap_bounderies((window_start, window_stop), read_interval)

                        coverage_low_bound = (window_length - 1) - abs(intersection[0] - window_start)
                        coverage_high_bound = (window_length - 1) - abs(intersection[1] - window_start)

                        read_length = read_interval[2]
                        if read_length not in coverage_mapping[chromosome]:
                            coverage_mapping[chromosome][read_length] = np.zeros(window_length, dtype=np.int)

                        if mapping_method == "fiveprime":
                            if read_interval[0] >= window_start:
                                coverage_mapping[chromosome][read_length][coverage_low_bound] += 1
                        elif mapping_method == "threeprime":
                            if read_interval[1] <= window_stop:
                                coverage_mapping[chromosome][read_length][coverage_high_bound] += 1

                        elif mapping_method == "centered":
                            center = round((read_interval[0] + read_interval[1]) / 2)
                            center_coverage = (window_length - 1) - abs(center - window_start)
                            if center >= window_start and center <= window_stop:
                                coverage_mapping[chromosome][read_length][center_coverage] += 1

                        elif mapping_method == "global":
                            coverage_mapping[chromosome][read_length][coverage_low_bound:coverage_high_bound+1] += 1

                else:
                    window_start = start - positions_out_ORF
                    window_stop = start + positions_in_ORF - 1
                    for read_interval in read_interval_dict[(chromosome, strand)].find((window_start, window_stop)):
                        intersection = misc.get_overlap_bounderies((window_start, window_stop), read_interval)

                        read_length = read_interval[2]
                        if read_length not in coverage_mapping[chromosome]:
                            coverage_mapping[chromosome][read_length] = np.zeros(window_length, dtype=np.int)

                        coverage_low_bound = intersection[0] - window_start
                        coverage_high_bound = intersection[1] - window_start

                        if mapping_method == "fiveprime":
                            if read_interval[1] <= window_stop:
                                coverage_mapping[chromosome][read_length][coverage_high_bound] += 1
                        elif mapping_method == "threeprime":
                            if read_interval[0] >= window_start:
                                coverage_mapping[chromosome][read_length][coverage_low_bound] += 1

                        elif mapping_method == "centered":
                            center = round((read_interval[0] + read_interval[1]) / 2)
                            center_coverage = center - window_start
                            if center >= window_start and center <= window_stop:
                                coverage_mapping[chromosome][read_length][center_coverage] += 1

                        elif mapping_method == "global":
                            coverage_mapping[chromosome][read_length][coverage_low_bound:coverage_high_bound+1] += 1

    return coverage_mapping
