#!/usr/bin/env python

"""
Script to transform the samples.tsv file into excel format
"""
import argparse
import pandas as pd

import excel_utils as eu

def convert_to_xlsx(samples_file, output_file):
    """
    Convert the .tsv file to xlsx

    :param samples_file: The samples.tsv file
    :param output_file: The output xlsx file
    """

    samples_df = pd.read_csv(samples_file, comment="#", sep="\t")

    # Show the file name rather than the whole path. str[-1] rather than str[1]:
    # the latter picked the second path component, so an absolute path became
    # "data" and a bare file name became NaN.
    for column in ("fastqFile", "fastqFile2"):
        if column in samples_df.columns:
            present = samples_df[column].notna()
            samples_df.loc[present, column] = (
                samples_df.loc[present, column].str.split("/").str[-1]
            )


    sheets = {"samples" : samples_df}

    eu.excel_writer(output_file, sheets, [])

def main():
    """
    Main function
    """

    parser = argparse.ArgumentParser(description='convert samples.tsv to xlsx')
    parser.add_argument("-i", "--samples", action="store", dest="samples", required=True, help= "The samples tsv file.")
    parser.add_argument("-o", "--xlsx", action="store", dest="output_file_path", required=True, help= "The output xlsx file.")
    args = parser.parse_args()

    convert_to_xlsx(args.samples, args.output_file_path)

if __name__ == '__main__':
    main()
