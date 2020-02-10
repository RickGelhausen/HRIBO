#!/usr/bin/env python
#This script requires numpy and pysam
#conda create -n mapping -c conda-forge -c bioconda pysam numpy
#Input bam files need to be sorted and indexed with samtools
#conda create -n samtools -c conda-forge -c bioconda samtools
#./mapping.py --bam_path RIBO-A-1_sorted.bam --wiggle_file_path /home/egg/current/SPP2002-data --no_of_aligned_reads 10000 --min_no_of_aligned_reads 10000
import numpy as np
import pysam
import argparse
from collections import defaultdict
from pprint import pprint
import os

def init_write_wig(file_handle,library_name, strand):
    tid = "%s_%s" % (library_name, strand)
    file_handle.write(("track type=wiggle_0 name=\"%s\"\n" % (tid)))


def write_wig(file_handle, seq_id, mappings, discard_zeros=True, factor=1.0):
    file_handle.write("variableStep chrom=%s span=1\n" % (seq_id))
    file_handle.write("\n".join(["%s %s" % (pos + 1, mapping * factor) for pos, mapping in filter(lambda pos_and_cov: pos_and_cov[1] != 0.0, enumerate(mappings))]) + "\n")

def add_aln_mapping(mappings, strand_swap, entry, increment, start, end, clip_length):
    if ((entry.is_reverse is False and entry.is_read2 is False) or
            (entry.is_reverse is True and entry.is_read2 is True)):
        if strand_swap:
           mappings["reverse"][start:end] -= increment
        else:
           mappings["forward"][start:end] += increment
    else:
        if strand_swap:
           mappings["forward"][start:end] += increment
        else:
           mappings["reverse"][start:end] -= increment

def add_first_base_mapping(mappings, strand_swap, entry, increment, start, end, clip_length):
    if ((entry.is_reverse is False and entry.is_read2 is False) or
            (entry.is_reverse is True and entry.is_read2 is True)):
        if strand_swap:
           mappings["reverse"][end-1] -= increment
        else:
           mappings["forward"][start] += increment
    else:
        if strand_swap:
           mappings["forward"][start] += increment
        else:
           mappings["reverse"][end-1] -= increment

def add_last_base_mapping(mappings, strand_swap, entry, increment, start, end, clip_length):
    if ((entry.is_reverse is False and entry.is_read2 is False) or
            (entry.is_reverse is True and entry.is_read2 is True)):
        if strand_swap:
           mappings["reverse"][start] -= increment
        else:
           mappings["forward"][end-1] += increment
    else:
        if strand_swap:
           mappings["forward"][end-1] += increment
        else:
           mappings["reverse"][start] -= increment


def add_centered_mapping(mappings, strand_swap, entry, increment, start, end, clip_length):
    center_start = start + clip_length
    center_end = end - clip_length
    center_length = float(center_end - center_start)
    if center_length < 1.0:
        return
    if ((entry.is_reverse is False and entry.is_read2 is False) or
            (entry.is_reverse is True and entry.is_read2 is True)):
        if strand_swap:
           mappings["reverse"][center_start:center_end] -= (increment / center_length)
        else:
           mappings["forward"][center_start:center_end] += (increment / center_length)
    else:
        if strand_swap:
           mappings["forward"][center_start:center_end] += (increment / center_length)
        else:
           mappings["reverse"][center_start:center_end] -= (increment / center_length)

def wig_file_path_setter(path,library_name,normtype,strand):
    resultpath=(path + normtype + "/" + library_name + "." + normtype + "." + strand + ".wig")
    return(resultpath)

def select_mapping_add_function(mapping_style):
    if mapping_style == "first_base_only":
        return add_first_base_mapping
    elif mapping_style == "last_base_only":
        return add_last_base_mapping
    elif mapping_style == "centered":
        return add_centered_mapping
    else:
        return add_aln_mapping

def compute_seqid_mapping(bam_path,read_count_splitting,mapping_add_function,strand_swap,clip_length):
    bam = pysam.Samfile(bam_path)
    results = []
    for ref_seq, length in zip(bam.references, bam.lengths):
        mappings = {}        
        for strand in ["forward", "reverse"]:
            mappings[strand] = np.array([0.0] * length)
        for entry in bam.fetch(ref_seq):
            number_of_hits = dict(entry.tags)["NH"]
            start = entry.pos
            end = entry.aend
            if read_count_splitting is True:
                increment = 1.0 / float(number_of_hits)
            else:
                increment = 1.0
            mapping_add_function(mappings, strand_swap,entry, increment, start, end, clip_length)
        results.append((ref_seq, mappings))
    return results

def compute_wig(bam_path,wig_file_path,library_name,min_no_of_aligned_reads,no_of_aligned_reads,read_count_splitting=True, mapping_style="centered", clip_length=11, non_strand_specific=False, strand_swap=False):
    mapping_add_function = select_mapping_add_function(mapping_style)
    seqid_mapping = compute_seqid_mapping(bam_path,read_count_splitting,mapping_add_function,strand_swap,clip_length)
    raw_fw_handle=open(wig_file_path_setter(wig_file_path,library_name,"raw","forward"), "w")
    init_write_wig(raw_fw_handle,library_name,"forward")
    raw_rw_handle=open(wig_file_path_setter(wig_file_path,library_name,"raw","reverse"), "w")
    init_write_wig(raw_rw_handle,library_name,"reverse")
    min_fw_handle=open(wig_file_path_setter(wig_file_path,library_name,"min","forward"), "w")
    init_write_wig(min_fw_handle,library_name,"forward")
    min_rw_handle=open(wig_file_path_setter(wig_file_path,library_name,"min","reverse"), "w")
    init_write_wig(min_rw_handle,library_name,"reverse")
    mil_fw_handle=open(wig_file_path_setter(wig_file_path,library_name,"mil","forward"), "w")
    init_write_wig(mil_fw_handle,library_name,"forward")
    mil_rw_handle=open(wig_file_path_setter(wig_file_path,library_name,"mil","reverse"), "w")
    init_write_wig(mil_fw_handle,library_name,"reverse")
    for (seqid, mappings) in seqid_mapping:
            write_wig(raw_fw_handle, seqid, mappings["forward"])
            write_wig(raw_rw_handle, seqid, mappings["reverse"])
            write_wig(min_fw_handle, seqid, mappings["forward"],factor=min_no_of_aligned_reads/no_of_aligned_reads)
            write_wig(min_rw_handle, seqid, mappings["reverse"],factor=min_no_of_aligned_reads/no_of_aligned_reads)
            write_wig(mil_fw_handle, seqid, mappings["forward"],factor=1000000/no_of_aligned_reads)
            write_wig(mil_rw_handle, seqid, mappings["reverse"],factor=1000000/no_of_aligned_reads)
    raw_fw_handle.close()
    raw_rw_handle.close()
    min_fw_handle.close()
    min_rw_handle.close()
    mil_fw_handle.close()
    mil_rw_handle.close()


def main():
    parser = argparse.ArgumentParser(description='Create single nucleotide mapping file')
    parser.add_argument("--bam_path", help='Bam path')
    parser.add_argument("--wiggle_file_path", help='File path to wig')
    parser.add_argument("--library_name", default="Library name", help='Library name to be displayed in the wig file')
    parser.add_argument("--mapping_style", default="centered", help='Mapping style: global, first_base_only, last_base_only, centered')
    parser.add_argument("--no_of_aligned_reads_file_path", help='File with total number of aligned reads')
    parser.add_argument("--clip_length", type=int, default=11, help='Clip length for centered mapping')
    parser.add_argument("--min_no_of_aligned_reads_file_path", help='File with minimal number of aligned reads for all used libraries')
    args = parser.parse_args()
    no_of_aligned_reads_file = open(args.no_of_aligned_reads_file_path,"r")
    no_of_aligned_reads = int(no_of_aligned_reads_file.read())
    min_no_of_aligned_reads_file = open(args.min_no_of_aligned_reads_file_path,"r")
    min_no_of_aligned_reads = int(min_no_of_aligned_reads_file.read())
    #mappings = {}
    compute_wig(args.bam_path,args.wiggle_file_path,args.library_name,min_no_of_aligned_reads,no_of_aligned_reads,True,args.mapping_style,11,False,False)
if __name__ == '__main__':
    main()
