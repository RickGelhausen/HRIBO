#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import pysam


def count_mapped_reads(args):
    """
    count the total number of reads for each input file
    return: None (output to file)
    """

    output_string_mapped = ""
    # run over all input files
    for bamfile in args.bamfiles:
        alignment_file = pysam.AlignmentFile(bamfile)
        total_mapped = {}
        for read in alignment_file.fetch():
            # get required attributes
            flag = int(read.flag) # flag: 0-forward-strand 4-unmapped 16-reverse-strand
            reference_name = read.reference_name

            if flag in [0,16]:

                if reference_name in total_mapped:
                    total_mapped[reference_name] += 1
                else:
                    total_mapped[reference_name] = 1

        wildcard = os.path.splitext(os.path.basename(bamfile))[0]

        for key, val in total_mapped.items():
            output_string_mapped += "%s\t%s\t%s\n" % (wildcard, key, str(val))

    with open(args.out_mapped, "w") as f:
        f.write(output_string_mapped)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='temporary script to calculate total mapped reads.')
    parser.add_argument("-b", "--bam", nargs="*", dest="bamfiles", required=True, help= "Read sequence files: (.bam/.sam)")
    parser.add_argument("-m", "--out_mapped", action="store", dest="out_mapped", required=True, help= "Output containg the average read length for each bam file.")
    args = parser.parse_args()

    count_mapped_reads(args)

if __name__ == '__main__':
    main()