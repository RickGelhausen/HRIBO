#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import pysam
import collections

def amount_overlap(r1, r2):
    return max(0, min(r1[1], r2[1]) - max(r1[0], r2[0]))

def get_read_information(args):
    """
    count the total number of reads for each input file
    return: None (output to file)
    """

    wildcards = []
    read_info_dict = {}
    # run over all input files
    for bamfile in args.bamfiles:
        alignment_file = pysam.AlignmentFile(bamfile)
        wildcard = os.path.splitext(os.path.basename(bamfile))[0]
        wildcards.append(wildcard)

        for read in alignment_file.fetch():
            # get required attributes
            flag = int(read.flag)
            reference_name = read.reference_name
            if flag == 0:
                strand = "+"
            elif flag == 16:
                strand = "-"

            if flag in [0,16]:
                start = read.pos + 1 # 0 based position
                length = len(read.query_sequence)
                stop = start + length

            if (wildcard, reference_name, strand) in read_info_dict:
                read_info_dict[(wildcard, reference_name, strand)].append((start, stop))
            else:
                read_info_dict[(wildcard, reference_name, strand)] = [(start, stop)]

    return read_info_dict, wildcards

def calculate_average_lengths(args, wildcards, read_info_dict):
    """
    for each row in the annotation, calculate the average mapped read length for each bam file
    """
    annotation_df = pd.read_csv(args.annotation, comment="#", sep="\t", header=None)

    column_count = 4 + len(wildcards)
    name_list = ["s%s" % str(x) for x in range(column_count)]
    nTuple = collections.namedtuple('Pandas', name_list)

    rows=[]
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        start = getattr(row, "_1")
        stop = getattr(row, "_2")
        strand = getattr(row, "_3")

        interval = (start, stop)
        average_lengths = []
        for wildcard in wildcards:
            read_count = 0
            total_overlap = 0
            for read in read_info_dict[(wildcard, reference_name, strand)]:
                overlap = amount_overlap(interval, read)
                if overlap > 0:
                    read_count += 1
                    total_overlap += overlap

            if read_count != 0:
                average_lengths.append(total_overlap/read_count)
            else:
                average_lengths.append(0)

        result = [reference_name, start, stop, strand] + average_lengths
        rows.append(nTuple(*result))

    out_df = pd.DataFrame.from_records(rows, columns=[x for x in range(column_count)])
    out_df.to_csv(args.out_length, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)



def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='temporary script to calculate total mapped reads.')
    parser.add_argument("-b", "--bam", nargs="*", dest="bamfiles", required=True, help= "Read sequence files: (.bam/.sam)")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "Annotation file to be filled with average_read_lengths.")
    parser.add_argument("-l", "--out_length", action="store", dest="out_length", required=True, help= "Output containg the average read length region and bam file.")
    args = parser.parse_args()

    read_info_dict, wildcards = get_read_information(args)
    calculate_average_lengths(args, wildcards, read_info_dict)


if __name__ == '__main__':
    main()
