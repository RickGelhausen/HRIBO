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

    output_string = ""
    # run over all input files
    for bamfile in args.bamfiles:
        alignment_file = pysam.AlignmentFile(bamfile)
        total_mapped = 0
        for read in alignment_file.fetch():
            # get required attributes
            flag = int(read.flag) # flag: 0-forward-strand 4-unmapped 16-reverse-strand

            if flag in [0,16]:
                total_length += len(read.query_sequence)
                total_mapped += 1

        wildcard = os.path.splitext(os.path.basename(bamfile))[0]
        output_string_mapped += "%s\t%s\n" % (wildcard, str(total_mapped))
        output_string_length += "%s\t%s\n" % (wildcard, str(total_length / total_mapped))

    with open(args.out_mapped, "w") as f:
        f.write(output_string_mapped)

    with open(args.out_length, "w") as f:
        f.write(output_string_length)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='temporary script to calculate total mapped reads.')
    parser.add_argument("-b", "--bam", nargs="*", dest="bamfiles", required=True, help= "Read sequence files: (.bam/.sam)")
    parser.add_argument("-m", "--out_mapped", action="store", dest="out_mapped", required=True, help= "Output containg the average read length for each bam file.")
    parser.add_argument("-l", "--out_length", action="store", dest="out_length", required=True, help= "Output containg total mapped reads for each bam file.")
    args = parser.parse_args()

    count_mapped_reads(args)

if __name__ == '__main__':
    main()
