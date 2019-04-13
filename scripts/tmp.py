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
                total_mapped += 1

        wildcard = os.path.splitext(os.path.basename(bamfile))[0]
        output_string += "%s\t%s\n" % (wildcard, str(total_mapped))

    with open(args.output, "w") as f:
        f.write(output_string)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='temporary script to calculate total mapped reads.')
    parser.add_argument("-b", "--bam", nargs="*", dest="bamfiles", required=True, help= "Read sequence files: (.bam/.sam)")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "Output containg total mapped reads for each bam file.")
    args = parser.parse_args()

    count_mapped_reads(args)

if __name__ == '__main__':
    main()
