#!/usr/bin/env python
import argparse
import os
import pysam
import csv
import matplotlib.pyplot as plt

def readlengthstats(input_bam_filepath,min_read_length,max_read_length,out_plot_filepath):
    readlengths = []
    readcounts = []
    readlength = min_read_length
    while readlength <= max_read_length:
        bamfile = pysam.AlignmentFile(input_bam_filepath, "rb")
        readcounter = 0
        for read in bamfile.fetch():
            if not (read.is_unmapped):
                if not (read.is_reverse):
                    currentreadlength = read.query_alignment_length
                    if currentreadlength == readlength:
                        readcounter += 1
        readlengths.append(readlength)
        readcounts.append(readcounter)
        print(str(readlength) + ":" + str(readcounter))
        readlength += 1
    print(readlengths)
    print(readcounts)
    plt.plot(readlengths, readcounts)
    plt.xlabel('read_lengths')
    plt.ylabel('counts')
    plt.savefig(out_plot_filepath + '/stat.pdf', format='pdf')


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='ReadLengthStat')
    parser.add_argument("--in_bam_filepath", help='Input bam path', required=True)
    parser.add_argument("--min_read_length", help='Minimal read length to consider', type=int, default=27)
    parser.add_argument("--max_read_length", help='Maximal read length to consider', type=int, default=33)
    parser.add_argument("--out_plot_filepath", help='Directory path to write output files, if not present the directory will be created', required=True)
    args = parser.parse_args()
    readlengthstats(args.in_bam_filepath,args.min_read_length,args.max_read_length,args.out_plot_filepath)


if __name__ == '__main__':
    main()
