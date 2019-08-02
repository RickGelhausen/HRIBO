#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections

def excel_writer(args, data_frames):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
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
            worksheet.set_column(idx, idx, max_len)
    writer.save()

def convert_to_xlsx(args):
    """
    covert the .tsv file to xlsx
    """

    nTuple = collections.namedtuple('Pandas', ["Method", "Condition", "Replicate", "Fastq file"])

    samples_df = pd.read_csv(args.samples, comment="#", sep="\t")
    samples_sheet = []
    for row in samples_df.itertuples(index=False, name='Pandas'):
        method = getattr(row, "method")
        condition = getattr(row, "condition")
        replicate = getattr(row, "replicate")
        fastq = getattr(row, "fastqFile")

        fastq = fastq.split("/")[1]
        samples_sheet.append(nTuple(method, condition, replicate, fastq))

    samples_df = pd.DataFrame.from_records(samples_sheet)

    sheets = {"samples" : samples_df}

    excel_writer(args, sheets)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='convert samples.tsv to xlsx')
    parser.add_argument("-i", "--samples", action="store", dest="samples", required=True, help= "the samples tsv file.")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    args = parser.parse_args()

    convert_to_xlsx(args)

if __name__ == '__main__':
    main()
