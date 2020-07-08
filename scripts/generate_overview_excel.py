#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections
import csv
from collections import Counter, OrderedDict
import itertools as iter

import interlap

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

import excel_utils as eu

class OrderedCounter(Counter, OrderedDict):
    pass

def create_interlap(annotation_dict):
    """
    create an interlap object to easily check for overlap of prediction and annotation
    """

    interlap_dict = {}

    tuples_dict = {}
    for key, val in annotation_dict.items():
        chrom, parts, strand = key.split(":")
        start, stop = parts.split("-")

        if (chrom, strand) in tuples_dict:
            tuples_dict[(chrom, strand)].append((int(start), int(stop), val[1]))
        else:
            tuples_dict[(chrom, strand)] = [(int(start), int(stop), val[1])]

    for key, val in tuples_dict.items():
        inter = interlap.InterLap()
        inter.update(val)
        interlap_dict[key] = inter

    return interlap_dict

def create_excel_file(args):
    # read annotation from file
    annotation_dict = eu.generate_annotation_dict(args.annotation_path)

    inter_dict = create_interlap(annotation_dict)

    # read the genome file
    genome_file = SeqIO.parse(args.genome_path, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    # get the total mapped reads for each bam file
    total_mapped_dict = {}
    with open(args.total_mapped, "r") as f:
        total = f.readlines()

    wildcards = []
    for line in total:
        wildcard, reference_name, value = line.strip().split("\t")
        total_mapped_dict[(wildcard, reference_name)] = int(value)
        wildcards.append(wildcard)

    wildcards = eu.get_unique(wildcards)

    TE_header = []
    for card in wildcards:
        if "RIBO" in card:
            TE_header.append("RIBO-"+card.split("-")[1])
        elif "TIS" in card and "RNATIS" not in card:
            TE_header.append("TIS-"+card.split("-")[1])

    counter = OrderedCounter(TE_header)

    TE_header = []
    for key, value in counter.items():
        for idx in range(value):
            TE_header.append("%s-%s" % (key,(idx+1)))
        if value > 1:
            TE_header.append("%s-avg" % key)

    conditions = []
    for card in wildcards:
        conditions.append(card.split("-")[1])

    conditions = eu.get_unique(conditions)

    xtail_dict, riborex_dict, deepribo_dict, reparation_dict = {}, {}, {}, {}
    if args.xtail_path != "":
        xtail_dict = eu.generate_xtail_dict(args.xtail_path)
    if args.riborex_path != "":
        riborex_dict = eu.generate_riborex_dict(args.riborex_path)
    if args.reads_deepribo != "":
        deepribo_dict = eu.generate_deepribo_dict(args.reads_deepribo)
    if args.reads_reparation != "":
        reparation_dict = eu.generate_reparation_dict(args.reads_reparation)

    #keys_union = list(set(deepribo_dict.keys()) | set(reparation_dict.keys()))

    keys_union = list(set().union(deepribo_dict.keys(), reparation_dict.keys(), annotation_dict.keys()))

    # read gff file
    all_sheet = []
    gff_file_rows = []

    contrasts = sorted(["%s-%s" %(tuple) for tuple in list(iter.combinations(conditions,2))])

    header = ["Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Overlapping_genes", "Old_locus_tag", "Name", "Gene_name", "Length", "Codon_count", "Start_codon", "Stop_codon"] +\
             ["Nucleotide_seq", "Aminoacid_seq"] +\
             [cond + "_TE" for cond in TE_header] +\
             [card + "_rpkm" for card in wildcards] +\
             ["Evidence_reparation", "Reparation_probability", "Evidence_deepribo", "Deepribo_rank", "Deepribo_score"] +\
             ["%s_%s" % (contrast, item) for contrast in contrasts for item in ["riborex_pvalue", "riborex_pvalue_adjusted", "riborex_log2FC"]] +\
             ["%s_%s" % (contrast, item) for contrast in contrasts for item in ["xtail_pvalue", "xtail_pvalue_adjusted", "xtail_log2FC"]]

    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)
    gffTuple = collections.namedtuple('Pandas', ["s1","s2","s3","s4","s5","s6","s7","s8","s9"])

    for key in keys_union:
        chromosome, mid, strand = key.split(":")
        start, stop = mid.split("-")

        locus_tag, locus_tag_overlap, name, gene_name = "","","",""
        reparation_probability, deepribo_rank, deepribo_score = 0, 0, 0
        evidence_reparation, evidence_deepribo = [], []

        read_list = []

        overlapping_intervals = list(inter_dict[(chromosome, strand)].find((int(start), int(stop))))

        if len(overlapping_intervals) > 0:
            locus_tag_overlap = ",".join([interval[2] for interval in overlapping_intervals])

        gene_id = ""
        old_locus_tag = ""
        if key in annotation_dict:
            gene_id, locus_tag, name, read_list, gene_name, old_locus_tag = annotation_dict[key]

        length = int(stop) - int(start) + 1
        codon_count = length / 3

        if key in reparation_dict:
            reparation_probability, reparation_evidence, read_list = reparation_dict[key]
            evidence_list = []
            for e in reparation_evidence.split(" "):
                if not "reparation" in e:
                    evidence_list.append("reparation-"+e)
                else:
                    evidence_list.append(e)
            evidence_reparation.extend(evidence_list)

        if key in deepribo_dict:
            deepribo_rank, deepribo_score, deepribo_evidence, read_list = deepribo_dict[key]
            evidence_list = []
            for e in deepribo_evidence.split(" "):
                if not "deepribo" in e:
                    evidence_list.append("deepribo-"+e)
                else:
                    evidence_list.append(e)
            evidence_deepribo.extend(evidence_list)

        xtail_list = []
        for contrast in contrasts:
            if (key,contrast) in xtail_dict:
                xtail_log2FC, xtail_pvalue, xtail_pvalue_adjusted = xtail_dict[(key,contrast)]
                xtail_list += [xtail_pvalue, xtail_pvalue_adjusted, xtail_log2FC]
            elif gene_id != "" and (gene_id,contrast) in xtail_dict:
                xtail_log2FC, xtail_pvalue, xtail_pvalue_adjusted = xtail_dict[(gene_id,contrast)]
                xtail_list += [xtail_pvalue, xtail_pvalue_adjusted, xtail_log2FC]
            else:
                xtail_list += [0,0,0]


        riborex_list = []
        for contrast in contrasts:
            if (key,contrast) in riborex_dict:
                riborex_log2FC, riborex_pvalue, riborex_pvalue_adjusted = riborex_dict[(key,contrast)]
                riborex_list += [riborex_pvalue, riborex_pvalue_adjusted, riborex_log2FC]
            elif gene_id != "" and (gene_id,contrast) in riborex_dict:
                riborex_log2FC, riborex_pvalue, riborex_pvalue_adjusted = riborex_dict[(gene_id,contrast)]
                riborex_list += [riborex_pvalue, riborex_pvalue_adjusted, riborex_log2FC]
            else:
                riborex_list += [0,0,0]

        start_codon, stop_codon, nucleotide_seq, aa_seq = eu.get_genome_information(genome_dict[chromosome], int(start)-1, int(stop)-1, strand)

        evidence_reparation.sort()
        evidence_deepribo.sort()

        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(eu.calculate_rpkm(total_mapped_dict[(wildcards[idx], chromosome)], val, length))

        TE_list = eu.calculate_TE(rpkm_list, wildcards, conditions)


        identifier = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
        evidence_reparation = " ".join(evidence_reparation)
        evidence_deepribo = " ".join(evidence_deepribo)
        if deepribo_rank == 0:
            deepribo_rank = "999999"
        result = [identifier, chromosome, start, stop, strand, locus_tag, locus_tag_overlap, old_locus_tag, name, gene_name, length, codon_count, start_codon, stop_codon] +\
                 [nucleotide_seq, aa_seq] +\
                 TE_list + rpkm_list +\
                 [evidence_reparation, reparation_probability, evidence_deepribo, deepribo_rank, deepribo_score]+\
                 riborex_list + xtail_list

        attributes = "ID=%s;" % identifier
        if locus_tag != "":
            attributes += "locus_tag=%s;" % locus_tag
        if old_locus_tag != "":
            attributes += "old_locus_tag=%s;" % old_locus_tag
        if gene_name != "":
            attributes += "Name=%s;" % gene_name
        elif name != "":
            attributes += "Name=%s;" % name
        else:
            attributes += "Name=%s;" % identifier

        gff_result = [chromosome, "HRIBO", "CDS", start, stop, ".", strand, ".", attributes]
        gff_file_rows.append(gffTuple(*gff_result))
        all_sheet.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(all_sheet, columns=[header[x] for x in range(len(header))])

    all_df = all_df.astype({"Start" : "int32", "Stop" : "int32"})
    all_df = all_df.sort_values(by=["Genome", "Start", "Stop"])

    all_df.to_csv(args.output_path.replace(".xlsx", ".tsv"), sep="\t", index=False, quoting=csv.QUOTE_NONE)

    dataframe_dict = { "all" : all_df }

    gff_df = pd.DataFrame.from_records(gff_file_rows, columns=[0,1,2,3,4,5,6,7,8])
    eu.excel_writer(args, dataframe_dict, wildcards)

    with open(args.output_path.replace(".xlsx", ".gff"), "w") as f:
        f.write("##gff-version 3\n")

    with open(args.output_path.replace(".xlsx", ".gff"), "a") as f:
        gff_df.to_csv(f, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create an overview table containing all results from the workflow.')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_path", required=True, help= "annotation file path.")
    parser.add_argument("-g", "--genome", action="store", dest="genome_path", required=True, help= "genome file path.")
    parser.add_argument("-r", "--riborex", action="store", dest="riborex_path", default="", help= "riborex csv file.")
    parser.add_argument("-x", "--xtail", action="store", dest="xtail_path", default="", help= "xtail csv file.")
    parser.add_argument("-o", "--xlsx", action="store", dest="output_path", required=True, help= "output xlsx file.")
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped", required=True, help= "file containing the total mapped reads for all alignment files.")
    parser.add_argument("--mapped_reads_deepribo", action="store", dest="reads_deepribo", default="", help= "file containing the individual read counts for deepribo.")
    parser.add_argument("--mapped_reads_reparation", action="store", dest="reads_reparation", default="", help= "file containing the individual read counts for reparation.")
    args = parser.parse_args()

    if args.reads_deepribo == "" and args.reads_reparation == "":
        sys.exit("no read files given!!")

    if args.xtail_path == "" and args.riborex_path == "":
        print("No differential expression data given. Proceeding without!")

    create_excel_file(args)

if __name__ == '__main__':
    main()
