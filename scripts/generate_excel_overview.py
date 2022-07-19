#!/usr/bin/env python
import argparse
import sys
import pandas as pd
import collections
import csv
import itertools as iter

import interlap

from Bio import SeqIO

import excel_utils as eu

FEATURE_MAP = { "ncrna" : "ncRNA", "trna" : "tRNA", "rrna" : "rRNA", "srna" : "sRNA",\
                "repeat_region" : "repeat_region", "start_codon" : "start_codons",\
                "stop_codon" : "stop_codon", "3'-utr" : "3'-UTR", "5'-utr" : "5'-UTR"}

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

def create_misc_excel_sheet(args, excel_sheet_dict, genome_dict, total_mapped_dict, wildcards, conditions, contrasts, te_header):
    """
    Create reduced sheets for all other features except CDS
    """
    annotation_dict =  eu.generate_non_cds_dict(args.annotation_path)

    xtail_dict = {}
    if args.xtail_path != "":
        xtail_dict = eu.generate_xtail_dict(args.xtail_path)
    riborex_dict = {}
    if args.riborex_path != "":
        riborex_dict = eu.generate_riborex_dict(args.riborex_path)
    deltate_dict = {}
    if args.deltate_path != "":
        deltate_dict = eu.generate_deltate_dict(args.deltate_path)

    sheet_row_dict = {}
    gff_rows = []

    header = ["Identifier", "Genome", "Start", "Stop", "Strand", "Feature", "Locus_tag", "Old_locus_tag", "Name", "Gene_name", "Length", "Codon_count", "Start_codon", "Stop_codon"] +\
             ["15nt upstream", "Nucleotide_seq", "Aminoacid_seq"] +\
             [f"{cond}_TE" for cond in te_header] +\
             [f"{card}_rpkm" for card in wildcards] +\
             [f"xtail_{contrast}_{item}" for contrast in contrasts for item in ["TE_log2FC", "TE_pvalue", "TE_pvalue_adjusted"]] +\
             [f"riborex_{contrast}_{item}" for contrast in contrasts for item in ["TE_log2FC", "TE_pvalue", "TE_pvalue_adjusted"]] +\
             [f"deltaTE_{contrast}_{item}" for contrast in contrasts for item in ["RIBO_log2FC", "RIBO_pvalue", "RIBO_pvalue_adjusted", "RNA_log2FC", "RNA_pvalue", "RNA_pvalue_adjusted", "TE_log2FC", "TE_pvalue", "TE_pvalue_adjusted"]] +\
             ["Product", "Note"]

    name_list = [f"s{x}" for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)
    gffTuple = collections.namedtuple('Pandas', ["s1","s2","s3","s4","s5","s6","s7","s8","s9"])

    for key in annotation_dict.keys():
        chromosome, mid, strand = key.split(":")
        start, stop = mid.split("-")

        feature, gene_id, locus_tag, name, gene_name, old_locus_tag, product, note, read_list = annotation_dict[key]

        length = int(stop) - int(start) + 1
        codon_count = length / 3

        xtail_list = []
        for contrast in contrasts:
            if (key, contrast) in xtail_dict:
                xtail_log2FC, xtail_pvalue, xtail_pvalue_adjusted = xtail_dict[(key, contrast)]
                xtail_list += [xtail_log2FC, xtail_pvalue, xtail_pvalue_adjusted]
            else:
                xtail_list += [None, None, None]

        riborex_list = []
        for contrast in contrasts:
            if (key, contrast) in riborex_dict:
                riborex_log2FC, riborex_pvalue, riborex_pvalue_adjusted = riborex_dict[(key, contrast)]
                riborex_list += [riborex_log2FC, riborex_pvalue, riborex_pvalue_adjusted]
            else:
                riborex_list += [None, None, None]

        deltate_list = []
        for contrast in contrasts:
            if (key, contrast) in deltate_dict:
                deltate_list += list(deltate_dict[(key, contrast)])
            else:
                deltate_list += [None, None, None, None, None, None, None, None, None]

        start_codon, stop_codon, nucleotide_seq, aa_seq, nt_window = eu.get_genome_information(genome_dict[chromosome], int(start)-1, int(stop)-1, strand)

        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(eu.calculate_rpkm(total_mapped_dict[(wildcards[idx], chromosome)], val, length))

        te_list = eu.calculate_te(rpkm_list, wildcards, conditions)

        identifier = f"{chromosome}:{start}-{stop}:{strand}"

        result = [identifier, chromosome, start, stop, strand, feature, locus_tag, old_locus_tag, name, gene_name, length, codon_count, start_codon, stop_codon] +\
                 [nt_window, nucleotide_seq, aa_seq] + te_list + rpkm_list + xtail_list + riborex_list + deltate_list + [product, note]

        attributes = f"ID={identifier};"
        if locus_tag != "":
            attributes += f"locus_tag={locus_tag};"
        if old_locus_tag != "":
            attributes += f"old_locus_tag={old_locus_tag};"
        elif name != "":
            attributes += f"Name={name};"
        else:
            attributes += f"Name={identifier};"

        gff_result = [chromosome, "HRIBO", feature, start, stop, ".", strand, ".", attributes]
        gff_rows.append(gffTuple(*gff_result))

        if feature.lower() in ["ncrna", "srna", "rrna", "trna", "start_codon", "stop_codon", "repeat_region", "3'-utr", "5'-utr"]:
            if FEATURE_MAP[feature.lower()] in sheet_row_dict:
                sheet_row_dict[FEATURE_MAP[feature.lower()]].append(nTuple(*result))
            else:
                sheet_row_dict[FEATURE_MAP[feature.lower()]] = [nTuple(*result)]
        else:
            if "misc" in sheet_row_dict:
                sheet_row_dict["misc"].append(nTuple(*result))
            else:
                sheet_row_dict["misc"] = [nTuple(*result)]

    gff_df = pd.DataFrame.from_records(gff_rows, columns=[0,1,2,3,4,5,6,7,8])

    with open(args.output_path.replace(".xlsx", "_misc.gff"), "w") as f:
        f.write("##gff-version 3\n")

    with open(args.output_path.replace(".xlsx", "_misc.gff"), "a") as f:
        gff_df.to_csv(f, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)

    for key in sheet_row_dict.keys():
        tmp_rows = sheet_row_dict[key]
        cur_df = pd.DataFrame.from_records(tmp_rows, columns=[header[x] for x in range(len(header))])

        cur_df = cur_df.astype({"Start" : "int32", "Stop" : "int32"})
        cur_df = cur_df.sort_values(by=["Genome", "Start", "Stop"])

        excel_sheet_dict[key] = cur_df

    return excel_sheet_dict

def create_cds_excel_sheet(args, excel_sheet_dict, genome_dict, total_mapped_dict, wildcards, conditions, contrasts, te_header):
    """
    Create the excel sheet for CDS.
    """

    # read annotation from file
    annotation_dict = eu.generate_annotation_dict(args.annotation_path)
    inter_dict = create_interlap(annotation_dict)

    xtail_dict, riborex_dict, deltate_dict, deepribo_dict, reparation_dict = {}, {}, {}, {}, {}
    if args.xtail_path != "":
        xtail_dict = eu.generate_xtail_dict(args.xtail_path)
    if args.riborex_path != "":
        riborex_dict = eu.generate_riborex_dict(args.riborex_path)
    if args.deltate_path != "":
        deltate_dict = eu.generate_deltate_dict(args.deltate_path)
    if args.reads_deepribo != "":
        deepribo_dict = eu.generate_deepribo_dict(args.reads_deepribo)
    if args.reads_reparation != "":
        reparation_dict = eu.generate_reparation_dict(args.reads_reparation)

    keys_union = list(set().union(deepribo_dict.keys(), reparation_dict.keys(), annotation_dict.keys()))

    # read gff file
    all_sheet = []
    annotated_sheet = []
    gff_file_rows = []

    header = ["Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Overlapping_genes", "Old_locus_tag", "Name", "Gene_name", "Length", "Codon_count", "Start_codon", "Stop_codon"] +\
             ["15nt upstream", "Nucleotide_seq", "Aminoacid_seq"] +\
             [f"{cond}_TE" for cond in te_header] +\
             [f"{card}_rpkm" for card in wildcards] +\
             ["Evidence_reparation", "Reparation_probability", "Evidence_deepribo", "Deepribo_rank", "Deepribo_score"] +\
             [f"xtail_{contrast}_{item}" for contrast in contrasts for item in ["TE_log2FC", "TE_pvalue", "TE_pvalue_adjusted"]] +\
             [f"riborex_{contrast}_{item}" for contrast in contrasts for item in ["TE_log2FC", "TE_pvalue", "TE_pvalue_adjusted"]] +\
             [f"deltaTE_{contrast}_{item}" for contrast in contrasts for item in ["RIBO_log2FC", "RIBO_pvalue", "RIBO_pvalue_adjusted", "RNA_log2FC", "RNA_pvalue", "RNA_pvalue_adjusted", "TE_log2FC", "TE_pvalue", "TE_pvalue_adjusted"]] +\
             ["Product", "Note"]

    name_list = [f"s{x}" for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)
    gffTuple = collections.namedtuple('Pandas', ["s1","s2","s3","s4","s5","s6","s7","s8","s9"])

    for key in keys_union:
        chromosome, mid, strand = key.split(":")
        start, stop = mid.split("-")

        locus_tag, locus_tag_overlap, name, gene_name, product, note = "","","","","",""
        reparation_probability, deepribo_rank, deepribo_score = 0, 0, 0
        evidence_reparation, evidence_deepribo = [], []

        read_list = []
        if (chromosome, strand) in inter_dict:
            overlapping_intervals = list(inter_dict[(chromosome, strand)].find((int(start), int(stop))))
        else:
            overlapping_intervals = []

        if len(overlapping_intervals) > 0:
            locus_tag_overlap = ",".join([interval[2] for interval in overlapping_intervals])

        gene_id = ""
        old_locus_tag = ""
        if key in annotation_dict:
            gene_id, locus_tag, name, gene_name, old_locus_tag, product, note, read_list = annotation_dict[key]

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
            if (key, contrast) in xtail_dict:
                xtail_log2FC, xtail_pvalue, xtail_pvalue_adjusted = xtail_dict[(key, contrast)]
                xtail_list += [xtail_log2FC, xtail_pvalue, xtail_pvalue_adjusted]
            else:
                xtail_list += [None, None, None]

        riborex_list = []
        for contrast in contrasts:
            if (key, contrast) in riborex_dict:
                riborex_log2FC, riborex_pvalue, riborex_pvalue_adjusted = riborex_dict[(key, contrast)]
                riborex_list += [riborex_log2FC, riborex_pvalue, riborex_pvalue_adjusted]
            else:
                riborex_list += [None, None, None]

        deltate_list = []
        for contrast in contrasts:
            if (key, contrast) in deltate_dict:
                deltate_list += list(deltate_dict[(key, contrast)])
            else:
                deltate_list += [None, None, None, None, None, None, None, None, None]

        start_codon, stop_codon, nucleotide_seq, aa_seq, nt_window = eu.get_genome_information(genome_dict[chromosome], int(start)-1, int(stop)-1, strand)

        evidence_reparation.sort()
        evidence_deepribo.sort()

        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(eu.calculate_rpkm(total_mapped_dict[(wildcards[idx], chromosome)], val, length))

        te_list = eu.calculate_te(rpkm_list, wildcards, conditions)

        identifier = f"{chromosome}:{start}-{stop}:{strand}"
        evidence_reparation = " ".join(evidence_reparation)
        evidence_deepribo = " ".join(evidence_deepribo)
        if deepribo_rank == 0:
            deepribo_rank = "999999"

        result = [identifier, chromosome, start, stop, strand, locus_tag, locus_tag_overlap, old_locus_tag, name, gene_name, length, codon_count, start_codon, stop_codon] +\
                 [nt_window, nucleotide_seq, aa_seq] +\
                 te_list + rpkm_list +\
                 [evidence_reparation, reparation_probability, evidence_deepribo, deepribo_rank, deepribo_score]+\
                  xtail_list + riborex_list + deltate_list + [product, note]

        attributes = f"ID={identifier};"
        if locus_tag != "":
            attributes += f"locus_tag={locus_tag};"
        if old_locus_tag != "":
            attributes += f"old_locus_tag={old_locus_tag};"
        if gene_name != "":
            attributes += f"Name={gene_name};"
        elif name != "":
            attributes += f"Name={name};"
        else:
            attributes += f"Name={identifier};"

        gff_result = [chromosome, "HRIBO", "CDS", start, stop, ".", strand, ".", attributes]
        gff_file_rows.append(gffTuple(*gff_result))
        all_sheet.append(nTuple(*result))
        if locus_tag != "" or name != "" or gene_name != "" or old_locus_tag != "":
            annotated_sheet.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(all_sheet, columns=[header[x] for x in range(len(header))])

    all_df = all_df.astype({"Start" : "int32", "Stop" : "int32"})
    all_df = all_df.sort_values(by=["Genome", "Start", "Stop"])

    all_df.to_csv(args.output_path.replace(".xlsx", ".tsv"), sep="\t", index=False, quoting=csv.QUOTE_NONE)

    annotated_df = pd.DataFrame.from_records(annotated_sheet, columns=[header[x] for x in range(len(header))])

    annotated_df = annotated_df.astype({"Start" : "int32", "Stop" : "int32"})
    annotated_df = annotated_df.sort_values(by=["Genome", "Start", "Stop"])

    excel_sheet_dict["all"] = all_df
    excel_sheet_dict["annotated"] = annotated_df

    gff_df = pd.DataFrame.from_records(gff_file_rows, columns=[0,1,2,3,4,5,6,7,8])

    with open(args.output_path.replace(".xlsx", ".gff"), "w") as f:
        f.write("##gff-version 3\n")

    with open(args.output_path.replace(".xlsx", ".gff"), "a") as f:
        gff_df.to_csv(f, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)

    return excel_sheet_dict

def create_excel_sheets(args):
    excel_sheet_dict = {}

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

    te_header = eu.get_te_header(wildcards)

    conditions = []
    for card in wildcards:
        conditions.append(card.split("-")[1])

    conditions = eu.get_unique(conditions)

    if args.contrasts is None or "," not in args.contrasts:
        contrasts = sorted([f"{x}-{y}" for x,y in list(iter.combinations(conditions, 2))], key= lambda s: s.lower())
    else:
        contrasts = args.contrasts.split(",")

    excel_sheet_dict = create_cds_excel_sheet(args, excel_sheet_dict, genome_dict, total_mapped_dict, wildcards, conditions, contrasts, te_header)
    excel_sheet_dict = create_misc_excel_sheet(args, excel_sheet_dict, genome_dict, total_mapped_dict, wildcards, conditions, contrasts, te_header)

    eu.excel_writer(args.output_path, excel_sheet_dict, wildcards)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create an overview table containing all results from the workflow.')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_path", required=True, help= "annotation file path.")
    parser.add_argument("-g", "--genome", action="store", dest="genome_path", required=True, help= "genome file path.")
    parser.add_argument("-x", "--xtail", action="store", dest="xtail_path", default="", help= "xtail csv file.")
    parser.add_argument("-r", "--riborex", action="store", dest="riborex_path", default="", help= "riborex csv file.")
    parser.add_argument("-d", "--deltate", action="store", dest="deltate_path", default="", help= "deltate csv file.")
    parser.add_argument("-o", "--xlsx", action="store", dest="output_path", required=True, help= "output xlsx file.")
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped", required=True, help= "file containing the total mapped reads for all alignment files.")
    parser.add_argument("-c", "--contrasts", action="store", dest="contrasts", default=None, help="file containing the contrasts for differential expression. If none provided, default sorting will be used.")
    parser.add_argument("--mapped_reads_deepribo", action="store", dest="reads_deepribo", default="", help= "file containing the individual read counts for deepribo.")
    parser.add_argument("--mapped_reads_reparation", action="store", dest="reads_reparation", default="", help= "file containing the individual read counts for reparation.")
    args = parser.parse_args()

    if args.reads_deepribo == "" and args.reads_reparation == "":
        sys.exit("no read files given!!")

    if args.xtail_path == "" and args.riborex_path == "" and args.deltate_path == "":
        print("No differential expression data given. Proceeding without!")

    create_excel_sheets(args)

if __name__ == '__main__':
    main()
