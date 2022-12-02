#!/usr/bin/env python
import pandas as pd
from pathlib import Path
import argparse
import os

def create_readcount_tables(read_count_file, contrast, tool, output_path):
    """
    Create a read count table for each method (RIBO and RNA) containing both conditions.
    """

    cond1, cond2 = contrast.split("-")

    if tool == "xtail":
        cond1, cond2 = cond2, cond1

    read_count_df = pd.read_csv(read_count_file, sep=",")

    ribo_cols = ["Identifier"] \
              + [col for col in read_count_df.columns if col.split("-")[:2] == ["RIBO", cond1]] \
              + [col for col in read_count_df.columns if col.split("-")[:2] == ["RIBO", cond2]]
    rna_cols = ["Identifier"] \
             + [col for col in read_count_df.columns if col.split("-")[:2] == ["RNA", cond1]] \
             + [col for col in read_count_df.columns if col.split("-")[:2] == ["RNA", cond2]]

    # create a read count table for each method
    ribo_df = read_count_df[ribo_cols]
    rna_df = read_count_df[rna_cols]

    ribo_df.to_csv(output_path / f"{contrast}_ribo_readcount_table.tsv", sep="\t", index=False)
    rna_df.to_csv(output_path / f"{contrast}_rna_readcount_table.tsv", sep="\t", index=False)

    # if tool == "riborex":
    #     condition_vector = [ "control" if col.split("-")[1] == cond1 else "treat" for col in ribo_cols[1:] ]
    # else:
    condition_vector = [ "treated" if col.split("-")[1] == cond1 else "control" for col in ribo_cols[1:] ]

    with open(output_path / f"{contrast}_condition_vector.csv", "w") as f:
        f.write(",".join(condition_vector) + "\n")

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Create the input files necessary to run riborex and xtail for a given contrast.')
    parser.add_argument("-r", "--read_count_file", action="store", dest="read_count_file", required=True, help="Read count table containing raw read counts.")
    parser.add_argument("-c", "--contrast", action="store", dest="contrast", required=True, help="Contrast to be used for preparaton.")
    parser.add_argument("-o", "--output_folder", action="store", dest="output_path", required=True, help= "The output folder.")
    parser.add_argument("-t", "--tool", action="store", dest="tool", required=True, help="Tool used (xtail/riborex).")
    args = parser.parse_args()
    p = Path(args.output_path)
    p.mkdir(parents=True, exist_ok=True)

    create_readcount_tables(args.read_count_file, args.contrast, args.tool, p)

if __name__ == '__main__':
    main()