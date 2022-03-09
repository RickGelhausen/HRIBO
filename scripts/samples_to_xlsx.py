#!/usr/bin/env python
import argparse
import pandas as pd

import excel_utils as eu

def convert_to_xlsx(args):
    """
    covert the .tsv file to xlsx
    """

    samples_df = pd.read_csv(args.samples, comment="#", sep="\t")

    column_names = [name.capitalize() for name in list(samples_df.columns)]

    samples_sheet = []
    for row in samples_df.itertuples(index=False, name='Pandas'):
        method = getattr(row, "method")
        condition = getattr(row, "condition")
        replicate = getattr(row, "replicate")
        fastq = getattr(row, "fastqFile")
        fastq = fastq.split("/")[1]
        ntup = [method, condition, replicate, fastq]

        if len(samples_df.columns) == 5:
            fastq2 = getattr(row, "fastqFile2")
            fastq2 = fastq2.split("/")[1]
            ntup.append(fastq2)

        samples_sheet.append(ntup)

    samples_df = pd.DataFrame.from_records(samples_sheet, columns=column_names)

    sheets = {"samples" : samples_df}

    eu.excel_writer(args.output, sheets, [])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='convert samples.tsv to xlsx')
    parser.add_argument("-i", "--samples", action="store", dest="samples", required=True, help= "the samples tsv file.")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    args = parser.parse_args()

    convert_to_xlsx(args)

if __name__ == '__main__':
    main()
