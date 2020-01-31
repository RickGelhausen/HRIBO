#!/usr/bin/env python
#This script requires numpy and pysam
#conda create -n coverage -c conda-forge -c bioconda pysam numpy
#Input bam files need to be sorted and indexed with samtools
#conda create -n samtools -c conda-forge -c bioconda samtools
#./coverage.py --bam_path RIBO-A-1_sorted.bam --wiggle_file_path /home/egg/current/SPP2002-data --no_of_aligned_reads 10000 --min_no_of_aligned_reads 10000
import numpy as np
import pysam
import argparse
from collections import defaultdict
     
def init_count_dict(ref_stats, refid):
    ref_stats[refid] = defaultdict(float)
    ref_stats[refid]["alns"]
    ref_stats[refid]["aligned_reads"]
    ref_stats[refid]["split_alns"]
    ref_stats[refid]["uniquely_aligned_reads"]
    ref_stats[refid]["aln_length_freqs"] = defaultdict(int)
    ref_stats[refid]["hits_read_freqs"] = defaultdict(int)

def calc_to_read(hits_read_freq):
    return dict((hits_read, freq/hits_read)
                for hits_read, freq in hits_read_freq.items())
        
def count_alignment(entry, refid, ref_stats,hits_read_freq):
    entry_tags_dict = dict(entry.tags)
    hits = entry_tags_dict["NH"]
    splits = float(entry_tags_dict.get("XL", 1))
    ref_stats[refid]["hits_read_freqs"][hits] += 1
    if "XL" in entry_tags_dict:
        ref_stats[refid]["split_alns"] += 1.0/splits
    ref_stats[refid]["alns"] += 1.0/splits
    ref_stats[refid]["aligned_reads"] += 1.0/(
        float(hits) * splits)
    if hits == 1:
        ref_stats[refid]["uniquely_aligned_reads"] += 1.0/splits
    ref_stats[refid]["aln_length_freqs"][entry.query_length] += 1

def sum_countings(ref_stats):
    total_stats = {}
    for refid, stats in ref_stats.items():
        for attribute, value in stats.items():
            if type(value) is int or type(value) is float:
                total_stats.setdefault(attribute, 0)
                total_stats[attribute] += value
            elif type(value) is dict:
                total_stats.setdefault(attribute, {})
                for value_int, freq in value.items():
                    total_stats[attribute].setdefault(value_int, 0)
                    total_stats[attribute][value_int] += freq
    return total_stats

def count_aln_reads_alns(stats, bam_path):
    bam = pysam.Samfile(bam_path)
    ref_stats = defaultdict(dict)
    hits_read_freq = {}
    for refid in bam.references:
        init_count_dict(ref_stats, refid)
    for entry in bam.fetch():
        refid = bam.get_reference_name(entry.tid)
        try:
            count_alignment(entry, refid, ref_stats, hits_read_freq)
        except KeyError:
            sys.stderr.write("Unspecified ref!\n")
            sys.exit(2)
    stats["stats_reference"] = ref_stats
    for refid, stats in ref_stats.items():
        ref_stats[refid]["hits_read_freqs"] = calc_to_read(ref_stats[refid]["hits_read_freqs"])
    stats["stats_total"] = sum_countings(ref_stats)
    return stats

        
def main():
    parser = argparse.ArgumentParser(description='Output read statistics')
    parser.add_argument("--bam_path", help='Bam path')
    args = parser.parse_args() 
    stats={}
    count_stats = count_aln_reads_alns(stats,args.bam_path)
    print(round(count_stats["stats_total"]["aligned_reads"]))
if __name__ == '__main__':
    main()
