#!/usr/bin/env python
# Input bam files need to be sorted and indexed with samtools
# conda create -n samtools -c conda-forge -c bioconda samtools
# ./mapping.py --bam_path RIBO-A-1_sorted.bam --wiggle_file_path /home/egg/current/SPP2002-data --no_of_aligned_reads 10000 --min_no_of_aligned_reads 10000

import argparse
import lib.library as hribo


def main():
    parser = argparse.ArgumentParser(description='Create single nucleotide mapping file')
    parser.add_argument("--bam_path", help='Bam path')
    parser.add_argument("--wiggle_file_path", help='File path to wig')
    parser.add_argument("--library_name", default="Library name", help='Library name to be displayed in the wig file')
    parser.add_argument("--mapping_style", default="centered", help='Mapping style: global, first_base_only, last_base_only, centered')
    parser.add_argument("--no_of_aligned_reads_file_path", help='File with total number of aligned reads')
    parser.add_argument("--clip_length", type=int, default=11, help='Clip length for centered mapping')
    args = parser.parse_args()
    #parse read count files
    (genome_read_dict,genome_min_read_dict) = hribo.get_read_count_dict(args.no_of_aligned_reads_file_path, args.library_name)
    #no_of_aligned_reads_file = open(args.no_of_aligned_reads_file_path,"r")
    #no_of_aligned_reads = int(no_of_aligned_reads_file.read())
    #min_no_of_aligned_reads_file = open(args.min_no_of_aligned_reads_file_path,"r")
    #min_no_of_aligned_reads = int(min_no_of_aligned_reads_file.read())
    #mappings = {}
    hribo.compute_wig(args.bam_path, args.wiggle_file_path, args.library_name, genome_read_dict, genome_min_read_dict, True, args.mapping_style, 11, False, False)


if __name__ == '__main__':
    main()
