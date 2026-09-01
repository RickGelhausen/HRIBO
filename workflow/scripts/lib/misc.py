"""
Contains miscellaneous scripts for different calculations.
Author: Rick Gelhausen
"""
import pandas as pd
import numpy as np

def get_overlap_bounderies(a, b):
    """
    get the overlap between two intervals
    """
    return (max(a[0], b[0]), min(a[1], b[1]))

def calculate_rpkm(gene_length, read_counts, total_counts):
    """
    Calculate rpkm for a gene
    """

    # rpm = read_counts / total_counts * 1e6
    # return rpm / gene_length * 1e3
    return (read_counts * 1000000000) / (total_counts * gene_length)

def count_reads(read_intervals_dict, chromosome, strand, beginning, end, mapping_method):
    """
    Select the right reads based on the mapping method.
    """
    read_intervals = read_intervals_dict[(chromosome,strand)].find((beginning, end))
    allowed_reads = []
    for read_interval in read_intervals:
        if mapping_method == "global":
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

    for chrom in metagene_coverage_dict:
        for read_length in metagene_coverage_dict[chrom]:
            metagene_coverage_dict[chrom][read_length] = (metagene_coverage_dict[chrom][read_length] / total_counts_dict[chrom]) * 1000000

    return metagene_coverage_dict

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

def create_data_frame(metagene_dict, positions_out_ORF, positions_in_ORF, state):
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
        dataframe_dict["no_evidence"] = pd.DataFrame({"coordinates": coordinates})
        return dataframe_dict

    for chrom in metagene_dict:
        if chrom not in dataframe_dict:
            dataframe_dict[chrom] = pd.DataFrame()
            dataframe_dict[chrom]["coordinates"] = coordinates
        for read_length in sorted(metagene_dict[chrom].keys()):
            dataframe_dict[chrom][f"{read_length}"] = metagene_dict[chrom][read_length]

    return dataframe_dict

def equalize_dictionary_keys(start_dict, stop_dict, positions_out_ORF, positions_in_ORF):
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

    read_lengths = {
        int(read_length)
        for coverage_dict in (start_dict, stop_dict)
        for chromosome_coverage in coverage_dict.values()
        for read_length in chromosome_coverage
    }
    if not read_lengths:
        return start_dict, stop_dict

    overall_min = min(read_lengths)
    overall_max = max(read_lengths)

    for key in unique_keys:
        for coverage_dict in (start_dict, stop_dict):
            for read_length in range(overall_min, overall_max + 1):
                if read_length not in coverage_dict[key]:
                    coverage_dict[key][read_length] = empty_window()

    return start_dict, stop_dict
