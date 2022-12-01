#!/usr/bin/env python
from calendar import c
import pandas as pd
from pathlib import Path
import argparse
import os

def split_bam_files(bam_folder, contrast):
    """
    Split the input bam_folder entries into RIBO and RNA files.
    """

    files = [entry for entry in Path(bam_folder).glob("*.bam") if entry.is_file()]
    files = sorted(files, key=lambda s: str(s.stem).lower())

    ribo_bam = []
    rna_bam = []
    contrasts = contrast.split("-")
    #contrasts = contrast[::-1]
    for contrast in contrasts:
        for file in files:
            file_prefix = file.stem
            method, condition, replicate = file_prefix.split("-")
            if condition != contrast:
                continue

            if "ribo" == method.lower():
                ribo_bam.append(file)

            if "rna" == method.lower():
                rna_bam.append(file)

    return ribo_bam, rna_bam

def prepare_read_count_dict(read_count_file):
    """
    Read all read counts from file.
    """

    read_count_df = pd.read_csv(read_count_file)

    return read_count_df.to_dict("series")

def compute_condition_map(ribo_bam, rna_bam):
    """
    Return a map of condition to integer
    """

    counter = 1
    condition_map = {}
    all_files = ribo_bam + rna_bam
    for file in all_files:
        _, condition, _ = file.stem.split("-")
        if condition not in condition_map:
            condition_map[condition] = counter
            counter += 1

    return condition_map

def create_sample_sheet(ribo_bam, rna_bam, output_path):
    """
    Create a sample sheet with RIBO and RNA samples
    """

    output_file = os.path.join(output_path, "samples_info.txt")
    replicate_file = os.path.join(output_path, "has_replicates.txt")

    condition_map = compute_condition_map(ribo_bam, rna_bam)

    replicate_dict = {}
    with open(output_file, "w") as of:
        of.write("SampleID\tCondition\tSeqType\tBatch\n")
        for ribo_file in ribo_bam:
            file_prefix = ribo_file.stem
            method, condition, replicate = file_prefix.split("-")
            if (method, condition) in replicate_dict:
                replicate_dict[(method, condition)].append(replicate)
            else:
                replicate_dict[(method, condition)] = [replicate]
            of.write(f"{file_prefix}\t{condition_map[condition]}\t{method}\t{replicate}\n")

        for rna_file in rna_bam:
            file_prefix = rna_file.stem
            method, condition, replicate = file_prefix.split("-")
            if (method, condition) in replicate_dict:
                replicate_dict[(method, condition)].append(replicate)
            else:
                replicate_dict[(method, condition)] = [replicate]
            of.write(f"{file_prefix}\t{condition_map[condition]}\t{method}\t{replicate}\n")

    has_replicates = True
    for key, val in replicate_dict.items():
        if len(val) < 2:
            has_replicates = False

    with open(replicate_file, "w") as of:
        of.write(f"{has_replicates}")

def create_readcount_table(bam_files, read_count_dict, output_path, file_name):
    """
    Create a read count file
    """

    output_file = os.path.join(output_path, file_name)

    columns = [file.stem for file in bam_files]

    counts_dict = {}
    counts_dict["Identifier"] = read_count_dict["Identifier"]
    for idx in range(len(columns)):
        counts_dict[columns[idx]] = read_count_dict[columns[idx]]

    read_counts_df = pd.DataFrame.from_dict(counts_dict)
    columns = read_counts_df.columns

    with open(output_file, "w") as f:
        f.write("\t".join([f"{column}" for column in columns[1:]]) + "\n")

    with open(output_file, "a") as f:
        read_counts_df.to_csv(f, sep="\t", header=False, index=False)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Create the input files necessary to run deltaTE for a given contrast.')
    parser.add_argument("-b", "--bam_folder", action="store", dest="bam_folder", required=True, help= "The folder containing all bam files.")
    parser.add_argument("-r", "--read_count_file", action="store", dest="read_count_file", required=True, help="Read count table containing raw read counts.")
    parser.add_argument("-c", "--contrast", action="store", dest="contrast", required=True, help="Contrast to be used for preparaton.")
    parser.add_argument("-o", "--output_folder", action="store", dest="output_path", required=True, help= "The output folder.")
    args = parser.parse_args()
    p = Path(args.output_path)
    p.mkdir(parents=True, exist_ok=True)

    ribo_bam, rna_bam = split_bam_files(args.bam_folder, args.contrast)
    create_sample_sheet(ribo_bam, rna_bam, args.output_path)

    read_count_dict = prepare_read_count_dict(args.read_count_file)
    create_readcount_table(ribo_bam, read_count_dict, args.output_path, "ribo_counts.txt")
    create_readcount_table(rna_bam, read_count_dict, args.output_path, "rna_counts.txt")

if __name__ == '__main__':
    main()