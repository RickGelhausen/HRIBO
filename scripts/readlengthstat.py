#!/usr/bin/env python
import argparse
import os
import pysam
import csv
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import ntpath
import json
import math
import lib.library as hribo


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='ReadLengthStat')
    parser.add_argument("--in_bam_filepath", help='Input bam path', required=True)
    parser.add_argument("--min_read_length", help='Minimal read length to consider', type=int, default=20)
    parser.add_argument("--max_read_length", help='Maximal read length to consider', type=int, default=50)
    parser.add_argument("--out_folder_filepath", help='Directory path to write output files, if not present the directory will be created', required=True)
    args = parser.parse_args()
    hribo.readlengthstats(args.in_bam_filepath,args.min_read_length,args.max_read_length,args.out_folder_filepath)


if __name__ == '__main__':
    main()
