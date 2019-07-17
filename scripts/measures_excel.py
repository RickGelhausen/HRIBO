#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

def retrieve_column_information(attributes):
    """
    check for gff2/gff3 format and generate a list of information for the final tables
    [locus_tag, name, product, note]
    """

    if ";" in attributes and "=" in attributes:
        attribute_list =  [x for x in re.split('[;=]', attributes)]
    else:
        attribute_list = [x.replace("\"", "") for x in re.split('[; ]', attributes) if x != ""]
 
    if len(attribute_list) % 2 == 0:
        for i in range(len(attribute_list)):
            if i % 2 == 0:
                attribute_list[i] = attribute_list[i].lower()    
    else:
        print(attributes)
        sys.exit("Attributes section of gtf/gff is wrongly formatted!")

    locus_tag = ""
    if "locus_tag" in attribute_list:
        locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

    name = ""
    if "name" in attribute_list:
        name = attribute_list[attribute_list.index("name")+1]
        if "cds" in name or "gene" in name:
            name = ""

    product = ""
    if "product" in attribute_list:
        product = attribute_list[attribute_list.index("product")+1]

    note = ""
    if "note" in attribute_list:
        note = attribute_list[attribute_list.index("note")+1]

    return [locus_tag, name, product, note]


def calculate_rpkm(total_mapped, read_count, read_length):
    """
    calculate the rpkm
    """
    return float("%.2f" % ((read_count * 1000000000) / (total_mapped * read_length)))

def calculate_tpm(normalize_factor, read_count, gene_length, average_length):
    """
    calculate the tpm
    """
    if read_count == 0:
        return 0

    return float("%.2f" % ((read_count * average_length * 1000000) / (normalize_factor * gene_length)))

def get_normalization_factor(read_df, wildcards, average_length_dict, prefix_columns):
    """
    calculate the denominator for the tpm
    """
    normalize_factor_dict = {}
    for row in read_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        start = int(getattr(row, "_1"))
        stop = int(getattr(row, "_2"))
        strand = getattr(row, "_5")

        gene_length = stop - start + 1

        read_list = [getattr(row, "_%s" %x) for x in range(prefix_columns,len(row))]
        for idx, val in enumerate(read_list):
            if (wildcards[idx], reference_name) in normalize_factor_dict:
                normalize_factor_dict[(wildcards[idx], reference_name)] += (int(val) * average_length_dict[(wildcards[idx],reference_name)]) / gene_length
            else:
                normalize_factor_dict[(wildcards[idx], reference_name)] = (int(val) * average_length_dict[(wildcards[idx],reference_name)]) / gene_length

    return normalize_factor_dict

def generate_excel_files(args):
    """
    generate and excel file for rpkm and one for tpm
    """
    # total_mapped_reads
    total_mapped_dict = {}
    with open(args.total_mapped, "r") as f:
        total = f.readlines()

    for line in total:
        wildcard, reference_name, value = line.strip().split("\t")
        total_mapped_dict[(wildcard, reference_name)] = int(value)

   # average read length
    average_length_dict = {}
    with open(args.length, "r") as f:
        total = f.readlines()

    for line in total:
        wildcard, reference_name, value = line.strip().split("\t")
        average_length_dict[(wildcard, reference_name)] = float(value)

    # read the comment containing the wildcards
    with open(args.reads, "r") as f:
        wildcards = f.readline()

    wildcards = wildcards.replace("# ", "").split()
    wildcards = [os.path.splitext(os.path.basename(card))[0] for card in wildcards]

    # read bed file
    read_df = pd.read_csv(args.reads, comment="#", header=None, sep="\t")

    rows_rpkm = []
    rows_tpm = []
    header_rpkm = ["Genome", "Source", "Feature", "Start", "Stop", "Strand", "Locus_tag", "Name", "Length", "Product", "Note"] + [card + "_rpkm" for card in wildcards]
    header_tpm = ["Genome", "Source", "Feature", "Start", "Stop", "Strand", "Locus_tag", "Name", "Length", "Product", "Note"] + [card + "_tpm" for card in wildcards]
    
    prefix_columns = len(read_df.columns) - len(wildcards)
    name_list = ["s%s" % str(x) for x in range(len(header_rpkm))]
    nTuple = collections.namedtuple('Pandas', name_list)
    normalize_factor_dict = get_normalization_factor(read_df, wildcards, average_length_dict, prefix_columns)
    for row in read_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        start = getattr(row, "_1")
        stop = getattr(row, "_2")
        strand = getattr(row, "_5")

        attributes = getattr(row, "_3")
        column_info = retrieve_column_information(attributes)

        source = getattr(row, "_6")
        feature = getattr(row, "_7")

        length = stop - start + 1
        read_list = [getattr(row, "_%s" %x) for x in range(prefix_columns,len(row))]

        ###################################### RPKM ##########################################
        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(calculate_rpkm(total_mapped_dict[(wildcards[idx], reference_name)], val, length))

        result_rpkm = [reference_name, source, feature, start, stop, strand, column_info[0], column_info[1], length, column_info[2], column_info[3]] + rpkm_list
        rows_rpkm.append(nTuple(*result_rpkm))

        ###################################### TPM ##########################################
        tpm_list = []
        for idx, val in enumerate(read_list):
            tpm_list.append(calculate_tpm(normalize_factor_dict[(wildcards[idx], reference_name)], val, length, average_length_dict[(wildcards[idx], reference_name)]))

        result_tpm = [reference_name, source, feature, start, stop, strand, column_info[0], column_info[1], length, column_info[2], column_info[3]] + tpm_list
        rows_tpm.append(nTuple(*result_tpm))

    excel_rpkm_df = pd.DataFrame.from_records(rows_rpkm, columns=[header_rpkm[x] for x in range(len(header_rpkm))])
    excel_tpm_df = pd.DataFrame.from_records(rows_tpm, columns=[header_tpm[x] for x in range(len(header_tpm))])

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
