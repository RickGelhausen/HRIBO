#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import collections


def excel_writer(args, data_frames, wildcards):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Note", "Aminoacid_seq", "Nucleotide_seq", "Start_codon", "Stop_codon", "Strand", "Codon_count"] + [card + "_rpkm" for card in wildcards]
    writer = pd.ExcelWriter(args.output, engine='xlsxwriter')
    for sheetname, df in data_frames.items():
        df.to_excel(writer, sheet_name=sheetname, index=False)
        worksheet = writer.sheets[sheetname]
        for idx, col in enumerate(df):
            series = df[col]
            if col in header_only:
                max_len = len(str(series.name)) + 2
            else:
                max_len = max(( series.astype(str).str.len().max(), len(str(series.name)) )) + 1
            print("Sheet: %s | col: %s | max_len: %s" % (sheetname, col, max_len))
            worksheet.set_column(idx, idx, max_len)
    writer.save()

def parse_orfs(args):
    # read the comment containing the wildcards
    with open(args.reads, "r") as f:
        wildcards = f.readline()

    wildcards = wildcards.replace("# ", "").split()
    wildcards = [os.path.splitext(os.path.basename(card))[0] for card in wildcards]

    #read bed file
    read_df = pd.read_csv(args.reads, comment="#", header=None, sep="\t")

    # read gff file
    main_sheet = []

    header = ["Class", "Feature count"] + [card + "_rpkm" for card in wildcards]
    prefix_columns = len(read_df.columns) - len(wildcards)
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    decode = {"srna" : "sRNA", "5'-utr" : "5'-UTR", "cds" : "CDS", "rrna" : "rRNA", "trna" : "tRNA", "transcript" : "transcript", "total" : "total"}
    feature_list = ["srna", "5'-utr", "cds", "rrna", "trna", "transcript", "total"]
    read_dict = collections.OrderedDict()
    count_dict = collections.OrderedDict()

    for f in feature_list:
        read_dict[f] = [0] * len(wildcards)
        count_dict[f] = 0

    main_sheet = []
    for row in read_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        start = getattr(row, "_1")
        stop = getattr(row, "_2")
        feature = getattr(row, "_7")

        read_list = [getattr(row, "_%s" %x) for x in range(prefix_columns,len(row))]
        if feature.lower() in read_dict:
            for idx, value in enumerate(read_list):
                print(value)
                read_dict[feature.lower()][idx] += value
                read_dict["total"][idx] += value
            count_dict[feature.lower()] += 1
            count_dict["total"] += 1
        else:
            print("feature not usable: " + feature)

    for key, val in read_dict.items():
        result = [decode[key], count_dict[key]] + val

        main_sheet.append(nTuple(*result))

    main_df = pd.DataFrame.from_records(main_sheet, columns=[header[x] for x in range(len(header))])

    dataframe_dict = { "Main" : main_df }

    excel_writer(args, dataframe_dict, wildcards)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel table with read information.')
    parser.add_argument("-r", "--mapped_reads", action="store", dest="reads", required=True, help= "file containing the individual read counts")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    args = parser.parse_args()

    parse_orfs(args)

if __name__ == '__main__':
    main()
