#!/usr/bin/env python
# samtools view -b -F 4 data/TIS-TolC-2.bam > data/TIS-TolC-2_mapped.bam
# samtools index data/TIS-TolC-2_mapped.bam
# samtools faidx exp5genome.fa
# time ./MetageneProfiling.py --in_bam_filepath ecoli/TIS-TolC-2.bam --in_gff_filepath ecoli/annotation.gff --out_plot_filepath ecoli/metagene
import argparse
import os
import numpy as np
import regex as re
import pysam
import csv
import matplotlib.pyplot as plt
import pandas as pd
import interlap
import random

from collections import defaultdict
from multiprocessing import Pool

def get_start_codons(input_gff_filepath):
    """
    read start codons from annotation file
    """
    startcodons = {"-": [], "+": []}
    seqid_set = set()
    with open(input_gff_filepath, newline='\n') as csvfile:
        gffreader = csv.reader(csvfile, delimiter='\t')
        for entry in gffreader:
            #filter empty and comment lines
            if (not len(entry) == 0 and not entry[0].startswith('#')):
                gfftype = entry[2]
                if gfftype == "CDS":
                    seqid = entry[0]
                    seqid_set.add(seqid)
                    source = entry[1]
                    gfftype = entry[2]
                    start = int(entry[3]) - 1
                    end = int(entry[4]) - 1
                    score = entry[5]
                    strand = entry[6]
                    phase = entry[7]
                    attributes = entry[8]
                    if strand == "+":
                        startcodons[strand].append(Feature("start", seqid, start, str(start+2), strand))
                    else:
                        startcodons[strand].append(Feature("start", seqid, str(end-2), end, strand))
    seqids = list(seqid_set)
    return seqids, startcodons


class Feature:
    def __init__(self, gfftype, seqid, start, end, strand):
        self.gfftype = gfftype
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand

    def __str__(self):
        return self.gfftype + "_" + self.seqid + "_" + self.start + "_" + self.end + "_" + self.strand

def get_overlap_bounderies(a, b):
    """
    get the overlap between two intervals
    """
    return (max(a[0], b[0]), min(a[1], b[1]))

def get_genome_sizes_dict(input_fai_filepath):
    """
    Read genome size from file and save them in a dictionary
    """
    genome_sizes_dict = {}
    with open(input_fai_filepath, newline='\n') as csvfile:
        tsvreader = csv.reader(csvfile, delimiter='\t')
        for entry in tsvreader:
            genome_sizes_dict[entry[0]] = int(entry[1])
    return(genome_sizes_dict)

#Metagene profiling
#For individual read lengths and indidividual strands and as summarized
#global, 5', 3', centered profiling for reads overlapping with area around the start codon

def meta_geneprofiling_p(in_gff_filepath, in_bam_filepath, out_plot_filepath, cpu_cores, min_read_length, max_read_length, normalization, in_fai_filepath):
    """
    Get all start_codons and all mapped reads,
    then run the metagene profiling on this data.
    """
    seqids, startcodons = get_start_codons(in_gff_filepath)
    genome_sizes_dict = get_genome_sizes_dict(in_fai_filepath)
    forward_length_reads_dict = defaultdict(list)
    reverse_length_reads_dict = defaultdict(list)
    bamfile = pysam.AlignmentFile(in_bam_filepath, "rb")
    if not os.path.exists(out_plot_filepath):
        os.makedirs(out_plot_filepath)
    print("Getting reads per length:")
    for read in bamfile.fetch():
        if not (read.is_unmapped):
            currentreadlength = read.query_length
            if currentreadlength >= min_read_length and currentreadlength <= max_read_length:
                if not read.is_reverse:
                    readfeature = Feature("read", read.reference_name, read.reference_start, read.reference_end, "+")
                    forward_length_reads_dict[currentreadlength].append(readfeature)
                else:
                    readfeature = Feature("read", read.reference_name, read.reference_start, read.reference_end, "-")
                    reverse_length_reads_dict[currentreadlength].append(readfeature)
    meta_gene_profiling(seqids, cpu_cores, forward_length_reads_dict, reverse_length_reads_dict, startcodons, out_plot_filepath, genome_sizes_dict, normalization)


def meta_gene_profiling(seqids, cpu_cores, forward_length_reads_dict, reverse_length_reads_dict, startcodons, out_plot_filepath, genome_sizes_dict, normalization):
    """
    Call the metagene mapping script for all seqids, read_lengths and mappings.
    Then plot all files.
    """
    for seqid in sorted(seqids):
        print("Metagene profiling for " + seqid + ":")
        chromosome_size = genome_sizes_dict.get(seqid, None)
        pool = Pool(processes=cpu_cores)
        threadsforward = pool.starmap(metagene_mapping, [(length, list(forward_length_reads_dict.get(length)), seqid, startcodons, "+") for length in sorted(forward_length_reads_dict)])
        globalforwardprofiles = dfObj = pd.DataFrame()
        fiveprimeforwardprofiles = dfObj = pd.DataFrame()
        threeprimeforwardprofiles = dfObj = pd.DataFrame()
        centeredforwardprofiles = dfObj = pd.DataFrame()
        if normalization:
                summarylabel = "mean"
        else:
                summarylabel = "sum"

        for (length, globalforwardmapping, fiveprimeforwardmapping, centeredforwardmapping, threeprimeforwardmapping) in threadsforward:
            globalforwardprofiles[length]=pd.Series(globalforwardmapping)
            fiveprimeforwardprofiles[length]=pd.Series(fiveprimeforwardmapping)
            centeredforwardprofiles[length]=pd.Series(centeredforwardmapping)
            threeprimeforwardprofiles[length]=pd.Series(threeprimeforwardmapping)

        globalforwardprofiles.loc[:,summarylabel] = globalforwardprofiles.sum(axis=1)
        fiveprimeforwardprofiles.loc[:,summarylabel] = fiveprimeforwardprofiles.sum(axis=1)
        centeredforwardprofiles.loc[:,summarylabel] = centeredforwardprofiles.sum(axis=1)
        threeprimeforwardprofiles.loc[:,summarylabel] = threeprimeforwardprofiles.sum(axis=1)
        plotprofile(globalforwardprofiles, seqid, out_plot_filepath, "global_forward", normalization)
        plotprofile(fiveprimeforwardprofiles, seqid, out_plot_filepath, "fiveprime_forward", normalization)
        plotprofile(centeredforwardprofiles, seqid, out_plot_filepath, "centered_forward", normalization)
        plotprofile(threeprimeforwardprofiles, seqid, out_plot_filepath, "threeprime_forward", normalization)
        threadsreverse = pool.starmap(metagene_mapping, [(length, list(reverse_length_reads_dict.get(length)), seqid, startcodons, "-") for length in sorted(reverse_length_reads_dict)])

        globalreverseprofiles = dfObj = pd.DataFrame()
        fiveprimereverseprofiles = dfObj = pd.DataFrame()
        threeprimereverseprofiles = dfObj = pd.DataFrame()
        centeredreverseprofiles = dfObj = pd.DataFrame()
        for (length, globalreversemapping, fiveprimereversemapping, centeredreversemapping, threeprimereversemapping) in threadsreverse:
            globalreverseprofiles[length]=pd.Series(globalreversemapping)
            fiveprimereverseprofiles[length]=pd.Series(fiveprimereversemapping)
            centeredreverseprofiles[length]=pd.Series(centeredreversemapping)
            threeprimereverseprofiles[length]=pd.Series(threeprimereversemapping)

        globalreverseprofiles.loc[:,summarylabel] = globalreverseprofiles.sum(axis=1)
        fiveprimereverseprofiles.loc[:,summarylabel] = fiveprimereverseprofiles.sum(axis=1)
        centeredreverseprofiles.loc[:,summarylabel] = centeredreverseprofiles.sum(axis=1)
        threeprimereverseprofiles.loc[:,summarylabel] = threeprimereverseprofiles.sum(axis=1)
        plotprofile(globalreverseprofiles, seqid, out_plot_filepath, "global_reverse", normalization)
        plotprofile(fiveprimereverseprofiles, seqid, out_plot_filepath, "fiveprime_reverse", normalization)
        plotprofile(centeredreverseprofiles, seqid, out_plot_filepath, "centered_reverse", normalization)
        plotprofile(threeprimereverseprofiles, seqid, out_plot_filepath, "threeprime_reverse", normalization)
        #merged plotting
        globalprofiles = globalforwardprofiles.add(globalreverseprofiles, fill_value=0)
        fiveprimeprofiles = fiveprimeforwardprofiles.add(fiveprimereverseprofiles, fill_value=0)
        centeredprofiles = centeredforwardprofiles.add(centeredreverseprofiles, fill_value=0)
        threeprimeprofiles = threeprimeforwardprofiles.add(threeprimereverseprofiles, fill_value=0)
        plotprofile(globalprofiles, seqid, out_plot_filepath, "global", normalization)
        plotprofile(fiveprimeprofiles, seqid, out_plot_filepath, "fiveprime", normalization)
        plotprofile(centeredprofiles, seqid, out_plot_filepath, "centered", normalization)
        plotprofile(threeprimeprofiles, seqid, out_plot_filepath, "threeprime", normalization)

def split_evenly(column_names, num_chunks):
    """
    split a list of column names in approximately equal sized chunks
    """
    k, m = divmod(len(column_names), num_chunks)
    return list((column_names[i * k + min(i, m):(i + 1) * k + min(i+1, m)] for i in range(num_chunks)))


def assign_color_list(data_split, color_list):
    custom_colors = []
    for ds in range(len(data_split)):
        custom_colors.append([color_list[x] for x in range(ds, len(color_list), len(data_split))])
    return custom_colors

def plotprofile(profiles, seqid, out_plot_filepath, profiletype, normalization):
    """
    Generate plot for the metagene profiling
    """

    if normalization == True:
        columns = list(profiles)
        for i in columns:
            total = profiles[i].sum()
            average_total = total / 500
            profiles[i] = profiles[i] / average_total

    profiles['coordinates'] = range(-100, -100 + len(profiles))
    profiles.to_csv(out_plot_filepath + "/" + seqid + "_" + profiletype + "_profiling.tsv", index=True, sep="\t", header=True,)
    profiles.to_excel(out_plot_filepath + "/" + seqid + "_" + profiletype + "_profiling.xlsx")

    column_names=list(profiles.columns)[:-2]

    #custom_colors = color_list[0:len(column_names[:-2])]
    cm = plt.get_cmap('gist_rainbow')
    color_list = [cm(1.*i/len(column_names)) for i in range(len(column_names))]
    max_Y = max(list(profiles.max(numeric_only=True))[:-2])

    if len(column_names) < 8:
        cur_ax = profiles.plot(x="coordinates", ylim=[0, max_Y + (max_Y * 5) / 100], color=color_list)
        cur_ax.set(xlabel="Position", ylabel="Coverage")
        cur_ax.axvline(x=0, color="grey")

        plt.savefig(out_plot_filepath + "/" + seqid + "_" + profiletype + "_profiling.pdf", format='pdf')
        plt.close()

    elif len(column_names) >= 8 and len(column_names) < 16:
        data_split = split_evenly(column_names, 2)

        custom_colors = assign_color_list(data_split, color_list)
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
        custom_axes = [ax1, ax2]
        for i in range(2):
            cur_ax = profiles.plot(ax=custom_axes[i], x="coordinates", y=data_split[i], ylim=[0, max_Y + (max_Y * 5) / 100], color=custom_colors[i])
            cur_ax.set(xlabel="Position", ylabel="Coverage")
            cur_ax.axvline(x=0, color="grey")

        plt.savefig(out_plot_filepath + "/" + seqid + "_" + profiletype + "_profiling.pdf", format='pdf')
        plt.close()

    else:
        data_split = split_evenly(column_names, 4)
        custom_colors = assign_color_list(data_split, color_list)
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))

        custom_axes = [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]
        for i in range(4):
            cur_ax = profiles.plot(ax=custom_axes[i], x="coordinates", y=data_split[i], ylim=[0, max_Y + (max_Y * 5) / 100], color=custom_colors[i])
            cur_ax.set(xlabel="Position", ylabel="Coverage")
            cur_ax.axvline(x=0, color="grey")

        plt.savefig(out_plot_filepath + "/" + seqid + "_" + profiletype + "_profiling.pdf", format='pdf')
        plt.close()

def metagene_mapping(length, length_reads, seqid, startcodons, strand):
    """
    For a given read length and seqid, calculate the coverage for all mapping methods
    """
    print("Len:" + str(length) + "\n")
    inter = interlap.InterLap()
    globalmapping = np.zeros(500, dtype=int)
    fiveprimemapping = np.zeros(500, dtype=int)
    centeredmapping = np.zeros(500, dtype=int)
    threeprimemapping = np.zeros(500, dtype=int)
    intervals = []
    for read in length_reads:
        if read.seqid == seqid:
            intervals.append((int(read.start), int(read.end)))


    if len(intervals) == 0:
        return length, globalmapping, fiveprimemapping, centeredmapping, threeprimemapping

    inter.update(intervals)
    for codon in startcodons[strand]:
        if strand == "+":
            codoninterval = (int(codon.start)-100, int(codon.end)+397)
        else:
            codoninterval = (int(codon.start)-397, int(codon.end)+100)
        overlapping_intervals = list(inter.find(codoninterval))
        if len(overlapping_intervals) == 0:
            continue
        for readinterval in overlapping_intervals:
            intersectiv = get_overlap_bounderies(codoninterval, readinterval)
            if strand =="+":
                coveragelower = intersectiv[0] - codoninterval[0]
                coverageupper = intersectiv[1] - codoninterval[0]
                if readinterval[0] >= codoninterval[0]:
                    fiveprimemapping[coveragelower] += 1
                if readinterval[1] <= codoninterval[1]:
                    threeprimemapping[coverageupper] += 1

                centerednumber = round((coveragelower + coverageupper)/2)
                centeredmapping[(centerednumber -1):(centerednumber +1)] += 1
                globalmapping[coveragelower:coverageupper] += 1

            else:
                coveragelower = 499 - abs(intersectiv[1] - codoninterval[0])
                coverageupper = 499 - abs(intersectiv[0] - codoninterval[0])
                if readinterval[1] <= codoninterval[1]:
                    fiveprimemapping[abs(coveragelower)] += 1
                if readinterval[0] >= codoninterval[0]:
                    threeprimemapping[abs(coverageupper)] += 1

                centerednumber = round((abs(coveragelower) + abs(coverageupper))/2)
                centeredmapping[(centerednumber -1):(centerednumber +1)] += 1
                globalmapping[abs(coveragelower):abs(coverageupper)] += 1

    return length, globalmapping, fiveprimemapping, centeredmapping, threeprimemapping

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='TISisitMetageneProfiling')
    parser.add_argument("--in_bam_filepath", help='Input bam path', required=True)
    parser.add_argument("--in_gff_filepath", help='Input gff path', required=True)
    parser.add_argument("--cpu_cores", help='Number of cpu cores to use', type=int, default=1)
    parser.add_argument("--min_read_length", help='Minimal read length to consider', type=int, default=27)
    parser.add_argument("--max_read_length", help='Maximal read length to consider', type=int, default=33)
    parser.add_argument("--out_plot_filepath", help='Directory path to write output files, if not present the directory will be created', required=True)
    parser.add_argument("--normalization", help='Toggles normalization by average read count per nucleotide', action='store_true')
    parser.add_argument("--in_fai_filepath", help='Input genome size (.fa.fai) path', required=True)
    args = parser.parse_args()
    meta_geneprofiling_p(args.in_gff_filepath, args.in_bam_filepath, args.out_plot_filepath, args.cpu_cores, args.min_read_length, args.max_read_length, args.normalization, args.in_fai_filepath)


if __name__ == '__main__':
    main()
