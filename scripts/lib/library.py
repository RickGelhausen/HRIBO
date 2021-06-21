import numpy as np
import pysam
import csv
import matplotlib.pyplot as plt
import pandas as pd
import interlap
from scipy.signal import find_peaks
from collections import defaultdict
from multiprocessing import Pool
import json
import statistics
import math
import os
import ntpath
from os import listdir
from os.path import isfile, join
import re
from pprint import pprint

def get_start_codons(input_gff_filepath):
    """
    read start codons from annotation file
    """
    startcodons = {"-": [], "+": []}
    seqid_set = set()
    with open(input_gff_filepath, newline='\n') as csvfile:
        gffreader = csv.reader(csvfile, delimiter='\t')
        for entry in gffreader:
            # filter empty and comment lines
            if (not len(entry) == 0 and not entry[0].startswith('#')):
                gfftype = entry[2]
                if gfftype == "CDS":
                    seqid = entry[0]
                    seqid_set.add(seqid)
                    # source = entry[1]
                    gfftype = entry[2]
                    start = int(entry[3]) - 1
                    end = int(entry[4]) - 1
                    # score = entry[5]
                    strand = entry[6]
                    # phase = entry[7]
                    # attributes = entry[8]
                    if strand == "+":
                        startcodons[strand].append(Feature("start", seqid, start, str(start+2), strand))
                    else:
                        startcodons[strand].append(Feature("start", seqid, str(end-2), end, strand))
    seqids = list(seqid_set)
    return seqids, startcodons


def get_stop_codons(input_gff_filepath):
    stopcodons = {"-": [], "+": []}
    seqid_set = set()
    with open(input_gff_filepath, newline='\n') as csvfile:
        gffreader = csv.reader(csvfile, delimiter='\t')
        for entry in gffreader:
            if (not entry[0].startswith('#')):
                gfftype = entry[2]
                if gfftype == "CDS":
                    seqid = entry[0]
                    seqid_set.add(seqid)
                    # source = entry[1]
                    gfftype = entry[2]
                    start = int(entry[3]) - 1
                    end = int(entry[4]) - 1
                    # score = entry[5]
                    strand = entry[6]
                    # phase = entry[7]
                    # attributes = entry[8]
                    if strand == "+":
                        stopcodons[strand].append(Feature("stop", seqid, str(end-2), end, strand))
                    else:
                        stopcodons[strand].append(Feature("stop", seqid, start, str(start+2), strand))
    seqids = list(seqid_set)
    return seqids, stopcodons


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

# Metagene profiling
# For individual read lengths and indidividual strands and as summarized
# global, 5', 3', centered profiling for reads overlapping with area around the start codon


def meta_geneprofiling_p(input_type, in_gff_filepath, in_bam_filepath, out_plot_filepath, cpu_cores, min_read_length, max_read_length, normalization, noise_reduction_analysis, in_readlengthstat_filepath):
    """
    Get all start_codons and all mapped reads,
    then run the metagene profiling on this data.
    """
    if input_type == "TIS":
        seqids, codons = get_start_codons(in_gff_filepath)
    else:
        seqids, codons = get_stop_codons(in_gff_filepath)
    if in_readlengthstat_filepath:
        if os.path.isfile(in_readlengthstat_filepath) and os.access(in_readlengthstat_filepath, os.R_OK):
            with open(in_readlengthstat_filepath) as json_file:
                length_reads_dict = json.load(json_file)
                max_reads_length = max(length_reads_dict, key=length_reads_dict.get)
                min_read_length = int(max_reads_length)-2
                max_read_length = int(max_reads_length)+2
    forward_length_reads_dict = defaultdict(list)
    reverse_length_reads_dict = defaultdict(list)
    bamfile = pysam.AlignmentFile(in_bam_filepath, "rb")
    if not os.path.exists(out_plot_filepath):
        os.makedirs(out_plot_filepath)
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
    meta_gene_profiling(input_type, seqids, cpu_cores, forward_length_reads_dict, reverse_length_reads_dict, codons, out_plot_filepath, normalization, noise_reduction_analysis)


def meta_gene_profiling(input_type, seqids, cpu_cores, forward_length_reads_dict, reverse_length_reads_dict, codons, out_plot_filepath, normalization, noise_reduction_analysis):
    """
    Call the metagene mapping script for all seqids, read_lengths and mappings.
    Then plot all files.
    """
    for seqid in sorted(seqids):
        pool = Pool(processes=cpu_cores)
        threadsforward = pool.starmap(metagene_mapping, [(length, list(forward_length_reads_dict.get(length)), seqid, codons, "+", input_type, noise_reduction_analysis) for length in sorted(forward_length_reads_dict)])
        globalforwardprofiles = dfObj = pd.DataFrame()
        fiveprimeforwardprofiles = dfObj = pd.DataFrame()
        threeprimeforwardprofiles = dfObj = pd.DataFrame()
        centeredforwardprofiles = dfObj = pd.DataFrame()
        if normalization:
            summarylabel = "mean"
        else:
            summarylabel = "sum"
        for (length, globalforwardmapping, fiveprimeforwardmapping, centeredforwardmapping, threeprimeforwardmapping) in threadsforward:
            globalforwardprofiles[length] = pd.Series(globalforwardmapping)
            fiveprimeforwardprofiles[length] = pd.Series(fiveprimeforwardmapping)
            centeredforwardprofiles[length] = pd.Series(centeredforwardmapping)
            threeprimeforwardprofiles[length] = pd.Series(threeprimeforwardmapping)
        threadsreverse = pool.starmap(metagene_mapping, [(length, list(reverse_length_reads_dict.get(length)), seqid, codons, "-", input_type, noise_reduction_analysis) for length in sorted(reverse_length_reads_dict)])
        globalreverseprofiles = dfObj = pd.DataFrame()
        fiveprimereverseprofiles = dfObj = pd.DataFrame()
        threeprimereverseprofiles = dfObj = pd.DataFrame()
        centeredreverseprofiles = dfObj = pd.DataFrame()
        for (length, globalreversemapping, fiveprimereversemapping, centeredreversemapping, threeprimereversemapping) in threadsreverse:
            globalreverseprofiles[length] = pd.Series(globalreversemapping)
            fiveprimereverseprofiles[length] = pd.Series(fiveprimereversemapping)
            centeredreverseprofiles[length] = pd.Series(centeredreversemapping)
            threeprimereverseprofiles[length] = pd.Series(threeprimereversemapping)
        globalforwardprofiles.loc[:,summarylabel] = globalforwardprofiles.sum(axis=1)
        fiveprimeforwardprofiles.loc[:,summarylabel] = fiveprimeforwardprofiles.sum(axis=1)
        centeredforwardprofiles.loc[:,summarylabel] = centeredforwardprofiles.sum(axis=1)
        threeprimeforwardprofiles.loc[:,summarylabel] = threeprimeforwardprofiles.sum(axis=1)
        globalreverseprofiles.loc[:,summarylabel] = globalreverseprofiles.sum(axis=1)
        fiveprimereverseprofiles.loc[:,summarylabel] = fiveprimereverseprofiles.sum(axis=1)
        centeredreverseprofiles.loc[:,summarylabel] = centeredreverseprofiles.sum(axis=1)
        threeprimereverseprofiles.loc[:,summarylabel] = threeprimereverseprofiles.sum(axis=1)
        #merged plotting
        globalprofiles = globalforwardprofiles.add(globalreverseprofiles, fill_value=0)
        fiveprimeprofiles = fiveprimeforwardprofiles.add(fiveprimereverseprofiles, fill_value=0)
        centeredprofiles = centeredforwardprofiles.add(centeredreverseprofiles, fill_value=0)
        threeprimeprofiles = threeprimeforwardprofiles.add(threeprimereverseprofiles, fill_value=0)
        plotprofile(globalforwardprofiles, seqid, out_plot_filepath, "global_forward", normalization, input_type, noise_reduction_analysis)
        plotprofile(fiveprimeforwardprofiles, seqid, out_plot_filepath, "fiveprime_forward", normalization, input_type, noise_reduction_analysis)
        plotprofile(centeredforwardprofiles, seqid, out_plot_filepath, "centered_forward", normalization, input_type, noise_reduction_analysis)
        plotprofile(threeprimeforwardprofiles, seqid, out_plot_filepath, "threeprime_forward", normalization, input_type, noise_reduction_analysis)
        plotprofile(globalreverseprofiles, seqid, out_plot_filepath, "global_reverse", normalization, input_type, noise_reduction_analysis)
        plotprofile(fiveprimereverseprofiles, seqid, out_plot_filepath, "fiveprime_reverse", normalization, input_type, noise_reduction_analysis)
        plotprofile(centeredreverseprofiles, seqid, out_plot_filepath, "centered_reverse", normalization, input_type, noise_reduction_analysis)
        plotprofile(threeprimereverseprofiles, seqid, out_plot_filepath, "threeprime_reverse", normalization, input_type, noise_reduction_analysis)
        plotprofile(globalprofiles, seqid, out_plot_filepath, "global", normalization, input_type, noise_reduction_analysis)
        plotprofile(fiveprimeprofiles, seqid, out_plot_filepath, "fiveprime", normalization, input_type, noise_reduction_analysis)
        plotprofile(centeredprofiles, seqid, out_plot_filepath, "centered", normalization, input_type, noise_reduction_analysis)
        plotprofile(threeprimeprofiles, seqid, out_plot_filepath, "threeprime", normalization, input_type, noise_reduction_analysis)


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


def find_optimal_offset(peaks,col):
    window_dict={}
    col_length=len(col)
    for peak in peaks:
        if (peak > 1) and (peak + 2 < col_length):
            window3_sum=math.floor(col[peak] + col[peak+1] + col[peak+2])
            if window3_sum>1:
                window_dict[window3_sum]=peak-100
        if peak > 2 and (peak + 1 < col_length):
            window2_sum=math.floor(col[peak-1] + col[peak] + col[peak+1])
            if window2_sum>1:
                window_dict[window2_sum]=peak-1-100
        if peak > 3 and (peak < col_length):
            window1_sum=math.floor(col[peak-2] + col[peak-1] + col[peak])
            if window1_sum>1:
                window_dict[window1_sum]=peak-2-100
    if window_dict:
        opt_window=max(window_dict.keys())
        return(window_dict.get(opt_window),opt_window)
    else:
        return("No peaks","no offset")


def plotprofile(profiles, seqid, out_plot_filepath, profiletype, normalization, input_type, noise_reduction_analysis):
    """
    Generate plot for the metagene profiling
    """
    (window_length, before_start_plus, after_start_plus, before_start_minus, after_start_minus) = set_window(input_type, noise_reduction_analysis)
    if normalization is True:
        columns = list(profiles)
        for i in columns:
            total = profiles[i].sum()
            average_total = total / window_length
            profiles[i] = profiles[i] / average_total
    profiles['coordinates'] = range(-before_start_plus, -before_start_plus + len(profiles))

    profiles.to_csv(out_plot_filepath + "/" + seqid + "_" + profiletype + "_profiling.tsv", index=True, sep="\t", header=True,)
    profiles.to_excel(out_plot_filepath + "/" + seqid + "_" + profiletype + "_profiling.xlsx")

    column_names=list(profiles.columns)[:-1]
    cm = plt.get_cmap('gist_rainbow')
    color_list = [cm(1.*i/len(column_names)) for i in range(len(column_names))]
    max_Y = max(list(profiles.max(numeric_only=True))[:-1])
    offset_dict = {}
    if len(column_names) < 8:
        # print(out_plot_filepath + "/" + seqid + "_" + profiletype)
        cur_ax = profiles.plot(x="coordinates", ylim=[0, max_Y + (max_Y * 5) / 100], color=color_list)
        cur_ax.set(xlabel="Position", ylabel="Coverage")
        cur_ax.axvline(x=0, color="grey")
        read_columns = column_names#[1:-1]# profiles.columns[1:-1]
        color_index=0
        # print(read_columns)
        average_read_length = math.floor(statistics.mean([x for x in read_columns if isinstance(x, (int,float))]))
        # average_read_length = math.floor(statistics.mean(filter(isInstance()columns if isinstance(x, int)]))
        for colname in read_columns:
            col =  profiles[colname]
            if (colname == "sum" or colname == "mean"):
                readlength = average_read_length
            else:
                readlength = int(colname)
            # peaks, properties=find_peaks(col.values, distance=readlength, width=[readlength-3,readlength+3])#,prominence=(0.1, 1))
            peaks, properties = find_peaks(col.values, distance=readlength, width=[1,40])#,prominence=(0.1, 1))
            peaks_offset =[math.floor(x) - 100 for x in peaks]
            xmins = [x - 100 for x in properties["left_ips"]]
            xmaxs = [x - 100 for x in  properties["right_ips"]]
            (optimal_offset,optimal_count) = find_optimal_offset(peaks,col)
            plt.hlines(y=properties["width_heights"], xmin=xmins, xmax=xmaxs, color=color_list[color_index])
            plt.vlines(x=peaks_offset, ymin=col[peaks] - properties["prominences"], ymax=col[peaks], color=color_list[color_index])
            color_index = color_index + 1
            offset_dict[colname] = str(optimal_offset) + "," + str(optimal_count)
        json_file_path = out_plot_filepath + "/" + seqid + "_" + profiletype + "_length_offset.json"
        json_file = open(json_file_path, 'w')
        json.dump(offset_dict, json_file)
        plt.savefig(out_plot_filepath + "/" + seqid + "_" + profiletype + "_profiling.pdf", format='pdf')
        plt.close()
    elif len(column_names) >= 8 and len(column_names) < 16:
        data_split = split_evenly(column_names, 2)
        custom_colors = assign_color_list(data_split, color_list)
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
        custom_axes = [ax1, ax2]
        if type(max_Y) == int or type(max_Y) == float:
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

def set_window(input_type,noise_reduction_analysis):
    #length of nucleotides to plot, before start + strand, after start + strand, before start - strand, after start minus strand
    if noise_reduction_analysis:
        d = { "TIS" : (200,100,97,97,100),
            "RIBO" : (200,100,97,97,100),
            "TTS" : (200,100,97,97,100)
            }
    else:
        d = { "TIS" : (500,100,397,397,100),
            "RIBO" : (500,100,397,397,100),
            "TTS" : (500,250,247,247,250)
            }
    window= d.get(input_type)
    if window == None:
        window=(500,100,397,397,100)
    return(window)

def metagene_mapping(length, length_reads, seqid, codons, strand, input_type, noise_reduction_analysis):
    """
    For a given read length and seqid, calculate the coverage for all mapping methods
    """
    inter = interlap.InterLap()
    (window_length,before_start_plus,after_start_plus,before_start_minus,after_start_minus) = set_window(input_type,noise_reduction_analysis)
    globalmapping = np.zeros(window_length, dtype=int)
    fiveprimemapping = np.zeros(window_length, dtype=int)
    centeredmapping = np.zeros(window_length, dtype=int)
    threeprimemapping = np.zeros(window_length, dtype=int)
    intervals = []
    for read in length_reads:
        if read.seqid == seqid:
            intervals.append((int(read.start), int(read.end)))

    if len(intervals) == 0:
        return length, globalmapping, fiveprimemapping, centeredmapping, threeprimemapping

    inter.update(intervals)
    for codon in codons[strand]:
        plus_strand_interval_begin = int(codon.start)-before_start_plus
        plus_strand_interval_end = int(codon.end)+after_start_plus
        minus_strand_interval_begin = int(codon.start)-before_start_minus
        minus_strand_interval_end = int(codon.end)+after_start_minus

        if strand == "+":
            codoninterval = (plus_strand_interval_begin, plus_strand_interval_end)
        else:
            codoninterval = (minus_strand_interval_begin, minus_strand_interval_end)

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
                coveragelower = (window_length -1) - abs(intersectiv[1] - codoninterval[0])
                coverageupper = (window_length -1) - abs(intersectiv[0] - codoninterval[0])
                if readinterval[1] <= codoninterval[1]:
                    fiveprimemapping[abs(coveragelower)] += 1
                if readinterval[0] >= codoninterval[0]:
                    threeprimemapping[abs(coverageupper)] += 1

                centerednumber = round((abs(coveragelower) + abs(coverageupper))/2)
                centeredmapping[(centerednumber -1):(centerednumber +1)] += 1
                globalmapping[abs(coveragelower):abs(coverageupper)] += 1

    return length, globalmapping, fiveprimemapping, centeredmapping, threeprimemapping


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
    peaks, properties=find_peaks(readcounts, width=[0,7])
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

def merge_offset(in_metagene_directorypath, out_path):
    # tis/tts
    offset_dict = {}
    profiling_type_dirs = [f.path for f in os.scandir(in_metagene_directorypath) if f.is_dir()]
    for profiling_type_dir in profiling_type_dirs:
        # norm/raw
        profiling_type_dir_key=os.path.basename(os.path.normpath(profiling_type_dir))
        offset_dict[profiling_type_dir_key]={}
        norm_type_dirs = [f.path for f in os.scandir(profiling_type_dir) if f.is_dir()]
        for norm_type_dir in norm_type_dirs:
            # samples
            norm_type_dir_key = os.path.basename(os.path.normpath(norm_type_dir))
            offset_dict[profiling_type_dir_key][norm_type_dir_key] = {}
            samples = [f.path for f in os.scandir(norm_type_dir) if f.is_dir()]
            for sample in samples:
                sample_key = os.path.basename(os.path.normpath(sample))
                offset_dict[profiling_type_dir_key][norm_type_dir_key][sample_key] = {}
                path = sample
                files = [f for f in listdir(sample) if isfile(join(sample, f))]
                json_files = [f for f in files if f.find(".json") != -1]
                json_wo_forward = [f for f in json_files if f.find("forward") == -1]
                json_wo_reverse_forward = [f for f in json_wo_forward if f.find("reverse") == -1]
                for json_file in json_wo_reverse_forward:
                    json_filename = os.path.basename(json_file)
                    json_filename = re.sub('\_length_offset.json$', '', json_filename)
                    offset_dict[profiling_type_dir_key][norm_type_dir_key][sample_key][json_filename] = {}
                    with open((path + "/" + json_file), 'r') as jf:
                        data = jf.read()
                        read_dict = json.loads(data)
                        offset_dict[profiling_type_dir_key][norm_type_dir_key][sample_key][json_filename] = read_dict
    # delete readlengthstat node obtained from adding length_reads_dict
    offset_dict.pop('readlengthstats', None)
    with open((out_path + "/" + 'merged_offsets.json'), 'w') as fp:
        json.dump(offset_dict, fp)
    return ""


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


def wig_file_path_setter(path, library_name, normtype, strand):
    resultpath = (path + normtype + "/" + library_name + "." + normtype + "." + strand + ".wig")
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


def compute_seqid_mapping(bam_path, read_count_splitting, mapping_add_function, strand_swap, clip_length):
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
            mapping_add_function(mappings, strand_swap, entry, increment, start, end, clip_length)
        results.append((ref_seq, mappings))
    return results


def compute_wig(bam_path, wig_file_path, library_name, genome_read_dict, genome_min_read_dict, read_count_splitting=True, mapping_style="centered", clip_length=11, non_strand_specific=False, strand_swap=False):
    mapping_add_function = select_mapping_add_function(mapping_style)
    seqid_mapping = compute_seqid_mapping(bam_path, read_count_splitting, mapping_add_function, strand_swap, clip_length)
    raw_fw_handle=open(wig_file_path_setter(wig_file_path, library_name, "raw", "forward"), "w")
    init_write_wig(raw_fw_handle,library_name,"forward")
    raw_rw_handle=open(wig_file_path_setter(wig_file_path, library_name, "raw", "reverse"), "w")
    init_write_wig(raw_rw_handle,library_name,"reverse")
    min_fw_handle=open(wig_file_path_setter(wig_file_path, library_name, "min", "forward"), "w")
    init_write_wig(min_fw_handle, library_name, "forward")
    min_rw_handle=open(wig_file_path_setter(wig_file_path, library_name, "min" , "reverse"), "w")
    init_write_wig(min_rw_handle,library_name,"reverse")
    mil_fw_handle=open(wig_file_path_setter(wig_file_path, library_name, "mil", "forward"), "w")
    init_write_wig(mil_fw_handle, library_name, "forward")
    mil_rw_handle=open(wig_file_path_setter(wig_file_path, library_name, "mil", "reverse"), "w")
    init_write_wig(mil_rw_handle, library_name, "reverse")
    for (seqid, mappings) in seqid_mapping:
        if seqid in genome_read_dict.keys():
            no_of_aligned_reads = int(genome_read_dict[seqid])
            min_no_of_aligned_reads = int(genome_min_read_dict[seqid])
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


def get_read_count_dict(input_read_filepath,library_name):
    genome_read_dict = {}
    genome_min_read_dict = {}
    with open(input_read_filepath, newline='\n') as csvfile:
        tsvreader = csv.reader(csvfile, delimiter='\t')
        for entry in tsvreader:
            if entry[0] == library_name:
               genome_read_dict[entry[1]] = int(entry[2])
            if entry[1] in genome_min_read_dict:
                if genome_min_read_dict[entry[1]] < entry[2]:
                    genome_min_read_dict[entry[1]] = entry[2]
            else:
                genome_min_read_dict[entry[1]] = entry[2]
    return(genome_read_dict, genome_min_read_dict)
