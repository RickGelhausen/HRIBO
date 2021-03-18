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

def readlengthstats(input_bam_filepath,min_read_length,max_read_length,out_folder_filepath):
    readlengths = []
    readcounts = []
    readlength = min_read_length
    length_count_dict={}
    count_length_dict={}
    while readlength <= max_read_length:
        bamfile = pysam.AlignmentFile(input_bam_filepath, "rb")
        readcounter = 0
        for read in bamfile.fetch():
            if not (read.is_unmapped):
                if not (read.is_reverse):
                    currentreadlength = read.query_alignment_length
                    if currentreadlength == readlength:
                        readcounter += 1
        length_count_dict[readlength]=readcounter
        count_length_dict[readcounter]=readlength
        readlengths.append(readlength)
        readcounts.append(readcounter)
        readlength += 1
    #print(length_count_dict)
    #print(count_length_dict)
    peaks, properties=find_peaks(readcounts,width=[1,7])
    #print(peaks)
    readlength_peaks=[readlengths[x] for x in peaks]
    #find max peak
    ymaxs=[length_count_dict.get(key) for key in readlength_peaks]
    index_max_peak=ymaxs.index(max(ymaxs))
    #print(readlength_peaks)
    plt.plot(readlengths, readcounts)
    xmins=[readlengths[int(math.floor(left))] for left in properties["left_ips"]]
    #print(xmins)
    xmaxs=[readlengths[int(math.ceil(right))] for right in properties["right_ips"]]
    #print(xmaxs)
    plt.hlines(y=properties["width_heights"], xmin = xmins, xmax = xmaxs)
    ymaxs=[length_count_dict.get(key) for key in readlength_peaks]
    #print(ymaxs)
    ymins =[0 for key in readlength_peaks]
    #print(ymins)
    max_length=readlength_peaks[index_max_peak]
    max_lower=xmins[index_max_peak]
    max_upper=xmaxs[index_max_peak]
    json_file_path=out_folder_filepath + ntpath.basename(input_bam_filepath) +  "_read_length_distribution.json"
    json_file = open(json_file_path, 'w')
    json.dump(length_count_dict, json_file)
    txt_file_path = out_folder_filepath + ntpath.basename(input_bam_filepath) +  "_read_length_distribution.txt"
    txt_file = open(txt_file_path, 'w')
    txt_file.write(str(max_length) + ":" + str(max_lower) + "-" + str(max_upper))
    txt_file.close()
    plt.vlines(x=readlength_peaks, ymin=ymins, ymax = ymaxs)
    plt.xlabel('read_lengths')
    plt.ylabel('counts')
    plt.savefig(out_folder_filepath + ntpath.basename(input_bam_filepath) +  "_read_length_distribution.pdf", format='pdf')


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='ReadLengthStat')
    parser.add_argument("--in_bam_filepath", help='Input bam path', required=True)
    parser.add_argument("--min_read_length", help='Minimal read length to consider', type=int, default=20)
    parser.add_argument("--max_read_length", help='Maximal read length to consider', type=int, default=50)
    parser.add_argument("--out_folder_filepath", help='Directory path to write output files, if not present the directory will be created', required=True)
    args = parser.parse_args()
    readlengthstats(args.in_bam_filepath,args.min_read_length,args.max_read_length,args.out_folder_filepath)


if __name__ == '__main__':
    main()
