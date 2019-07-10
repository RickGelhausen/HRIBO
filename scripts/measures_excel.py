#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import collections

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

def calculate_rpkm(total_mapped, read_count, read_length):
    """
    calculate the rpkm
    """
    return "%.2f" % ((read_count * 1000000000) / (total_mapped * read_length))

def calculate_tpm(normalize_factor, read_count, gene_length, average_length):
    """
    calculate the tpm
    """
    return "%.2f" % ((read_count * average_length * 1000000) / (normalize_factor * gene_length))

def get_normalization_factor(read_df, wildcards, average_length_dict):
    normalize_factor_dict = {}
    for wildcard in wildcards:
        normalize_factor_dict[wildcard] = 0

    for row in read_df.itertuples(index=False, name='Pandas'):
        start = int(getattr(row, "_1"))
        stop = int(getattr(row, "_2"))

        read_length = stop - start + 1

        read_list = [getattr(row, "_%s" %x) for x in range(5,len(row))]
        for idx, val in enumerate(read_list):
            normalize_factor_dict[wildcards[idx]] += (int(val) * average_length_dict[wildcards[idx]]) / read_length

    return normalize_factor_dict


def generate_excel_files(args):
    # total_mapped_reads
    total_mapped_dict = {}
    with open(args.total_mapped, "r") as f:
        total = f.readlines()

    for line in total:
        key, value = line.strip().split("\t")
        total_mapped_dict[key] = int(value)

    # average read length
    average_length_dict = {}
    with open(args.length, "r") as f:
        total = f.readlines()

    for line in total:
        key, value = line.strip().split("\t")
        average_length_dict[key] = float(value)

    # read the comment containing the wildcards
    with open(args.reads, "r") as f:
        wildcards = f.readline()

    wildcards = wildcards.replace("# ", "").split()
    wildcards = [os.path.splitext(os.path.basename(card))[0] for card in wildcards]

    read_df = pd.read_csv(args.reads, comment="#", header=None, sep="\t")
    column_count = len(read_df.columns)

    name_list = ["s%s" % str(x) for x in range(column_count)]
    nTuple = collections.namedtuple('Pandas', name_list)
    # read gff file
    rows_rpkm = []
    rows_tpm = []
    header = ["orfID", "start", "stop", "strand", "length"] + [card + "_rpkm" for card in wildcards]
    rows_rpkm.append(nTuple(*header))
    header = ["orfID", "start", "stop", "strand", "length"] + [card + "_tpm" for card in wildcards]
    rows_tpm.append(nTuple(*header))

    normalize_factor_dict = get_normalization_factor(read_df, wildcards, average_length_dict)
    for row in read_df.itertuples(index=False, name='Pandas'):
        id = getattr(row, "_0")
        start = getattr(row, "_1")
        stop = getattr(row, "_2")
        strand = getattr(row, "_3")

        length = stop - start + 1
        read_list = [getattr(row, "_%s" %x) for x in range(5,len(row))]

        ###################################### RPKM ##########################################
        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(calculate_rpkm(total_mapped_dict[wildcards[idx]], val, length))

        result_rpkm = [id, start, stop, strand, length] + rpkm_list
        rows_rpkm.append(nTuple(*result_rpkm))

        ###################################### TPM ##########################################
        tpm_list = []
        for idx, val in enumerate(read_list):
            tpm_list.append(calculate_tpm(normalize_factor_dict[wildcards[idx]], val, length, average_length_dict[wildcards[idx]]))

        result_tpm = [id, start, stop, strand, length] + tpm_list
        rows_tpm.append(nTuple(*result_tpm))

    excel_rpkm_df = pd.DataFrame.from_records(rows, columns=[x for x in range(column_count)])
    excel_tpm_df = pd.DataFrame.from_records(rows, columns=[x for x in range(column_count)])

    excel_rpkm_df.to_excel(args.rpkm, sheet_name=args.sheet_name)
    excel_tpm_df.to_excel(args.tpm, sheet_name=args.sheet_name)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel files containing: rpkm and tpm measures')
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped", required=True\
                                                    , help= "file containing the total mapped reads for all alignment files.")
    parser.add_argument("-r", "--mapped_reads", action="store", dest="reads", required=True, help= "file containing the individual read counts.")
    parser.add_argument("-l", "--read_lengths", action="store", dest="length", required=True, help= "file containing the average read lengths.")
    parser.add_argument("--out_rpkm", action="store", dest="rpkm", required=True, help= "output rpkm xlsx file")
    parser.add_argument("--out_tpm", action="store", dest="tpm", required=True, help= "output tpm xlsx file")
    parser.add_argument("-s", "--sheet_name", action="store", dest="sheet_name", help= "name of the excel sheet", default="Sheet 1")

    args = parser.parse_args()

    generate_excel_files(args)

if __name__ == '__main__':
    main()
