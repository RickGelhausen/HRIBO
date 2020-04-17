#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections
import csv
from collections import Counter, OrderedDict

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

class OrderedCounter(Counter, OrderedDict):
    pass

def calculate_rpkm(total_mapped, read_count, read_length):
    """
    calculate the rpkm
    """
    return float("%.2f" % ((read_count * 1000000000) / (total_mapped * read_length)))

def get_unique(in_list):
    seen = set()
    seen_add = seen.add
    return [x for x in in_list if not (x in seen or seen_add(x))]

def retrieve_column_information(attributes):
    """
    check for gff2/gff3 format and generate a list of information for the final tables
    [locus_tag, name, product, note]
    """

    if ";" in attributes and "=" in attributes:
        attribute_list = [x for x in re.split('[;=]', attributes) if x != ""]
    else:
        attribute_list = [x.replace(";", "") for x in list(csv.reader([attributes], delimiter=' ', quotechar='"'))[0]]

    if "ORF_type=;" in attributes:
        attribute_list.remove("ORF_type")

    if len(attribute_list) % 2 == 0:
        for i in range(len(attribute_list)):
            if i % 2 == 0:
                attribute_list[i] = attribute_list[i].lower()
    else:
        print(attributes)
        sys.exit("Attributes section of gtf/gff is wrongly formatted!")

    locus_tag = ""
    if "locus_tag" in attribute_list:
        locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

    name = ""
    if "name" in attribute_list:
        name = attribute_list[attribute_list.index("name")+1]

    return [locus_tag, name]

def get_genome_information(genome, start, stop, strand):
    """
    retrieve the nucleotide sequence and amino acid sequence
    and the start and stop codons
    """
    if strand == "+":
        nucleotide_seq = genome[0][start:stop+1]
    else:
        nucleotide_seq = genome[1][start:stop+1][::-1]

    start_codon = nucleotide_seq[0:3]
    stop_codon = nucleotide_seq[-3:]

    coding_dna = Seq(nucleotide_seq, generic_dna)
    if len(coding_dna) % 3 != 0:
        aa_seq = ""
    else:
        aa_seq = str(coding_dna.translate(table=11,to_stop=True))

    return start_codon, stop_codon, nucleotide_seq, aa_seq

def excel_writer(args, data_frames, wildcards):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Note", "Aminoacid_seq", "Nucleotide_seq", "Start_codon", "Stop_codon", "Strand", "Codon_count"] + [card + "_rpkm" for card in wildcards]
    writer = pd.ExcelWriter(args.output_path, engine='xlsxwriter')
    for sheetname, df in data_frames.items():
        df.to_excel(writer, sheet_name=sheetname, index=False)
        worksheet = writer.sheets[sheetname]
        for idx, col in enumerate(df):
            series = df[col]
            if col in header_only:
                max_len = len(str(series.name)) + 2
            else:
                max_len = max(( series.astype(str).str.len().max(), len(str(series.name)) )) + 1
            print("Sheet: %s | col: %s | max_len: %s" % (sheetname, col, max_len))
            worksheet.set_column(idx, idx, max_len)
    writer.save()

def TE(ribo_count, rna_count):
    """
    calculate the translational efficiency for one entry
    """

    if ribo_count == 0 and rna_count == 0:
        return 0
    elif rna_count == 0:
        return ribo_count
    else:
        return ribo_count / rna_count

def calculate_TE(read_list, wildcards, conditions):
    """
    calculate the translational efficiency
    """
    read_dict = OrderedDict()
    for idx in range(len(wildcards)):
        method, condition, replicate = wildcards[idx].split("-")
        key = (method, condition)
        if key in read_dict:
            read_dict[key].append(read_list[idx])
        else:
            read_dict[key] = [read_list[idx]]

    TE_list = []
    for cond in conditions:
        if ("RIBO", cond) in read_dict:
            if len(read_dict[("RIBO", cond)]) == len(read_dict[("RNA", cond)]):
                ribo_list = read_dict[("RIBO", cond)]
                rna_list = read_dict[("RNA", cond)]
                t_eff = [TE(ribo_list[idx],rna_list[idx]) for idx in range(len(ribo_list))]

                if len(t_eff) > 1:
                    t_eff.extend([sum(t_eff) / len(ribo_list)])

                t_eff = [float("%.2f" % x) for x in t_eff]
                TE_list.extend(t_eff)
            else:
                TE_list.extend([0])

        if ("TIS", cond) in read_dict:
            if len(read_dict[("TIS", cond)]) == len(read_dict[("RNATIS", cond)]):
                ribo_list = read_dict[("TIS", cond)]
                rna_list = read_dict[("RNATIS", cond)]

                t_eff = [TE(ribo_list[idx],rna_list[idx]) for idx in range(len(ribo_list))]

                if len(t_eff) > 1:
                    t_eff.extend([sum(t_eff) / len(ribo_list)])

                t_eff = [float("%.2f" % x) for x in t_eff]
                TE_list.extend(t_eff)
            else:
                TE_list.extend([0])

    return TE_list

def generate_riborex_dict(riborex_path):
    """
    create a dictionary containing all important riborex input
    {geneID : (log2fc, pvalue, pvalue_adj)}
    """

    riborex_df = pd.read_csv(riborex_path, sep=",", comment="#")
    riborex_dict = {}

    for row in riborex_df.itertuples(index=False, name='Pandas'):
        gene_id = getattr(row, "gene_id")
        log2fc = getattr(row, "log2FoldChange")
        pvalue = getattr(row, "pvalue")
        pvalue_adj = getattr(row, "padj")
        contrasts = getattr(row, "contrast")

        riborex_dict[gene_id] = (log2fc, pvalue, pvalue_adj, contrasts)

    return riborex_dict


def generate_xtail_dict(xtail_path):
    """
    create a dictionary containing all important xtail input
    {geneID : (log2fc, pvalue, pvalue_adj)}
    """

    xtail_df = pd.read_csv(xtail_path, sep=",", comment="#")
    xtail_dict = {}

    for row in xtail_df.itertuples(index=False, name='Pandas'):
        gene_id = getattr(row, "gene_id")
        log2fc = getattr(row, "log2FC_TE_final")
        pvalue = getattr(row, "pvalue_final")
        pvalue_adj = getattr(row, "pvalue_adjust")
        contrasts = getattr(row, "contrast")

        xtail_dict[gene_id] = (log2fc, pvalue, pvalue_adj, contrasts)

    return xtail_dict


def generate_reparation_dict(reparation_path):
    """
    create a dictionary containing all important reparation input
    """

    reparation_df = pd.read_csv(reparation_path, sep="\t", comment="#")
    prefix_columns = 9

    reparation_dict = {}
    for row in reparation_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attribute_list = [x for x in re.split('[;=]', getattr(row, "_8")) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")

        proba = ""
        if "prob" in attribute_list:
            proba = attribute_list[attribute_list.index("prob")+1]

        evidence = ""
        if "evidence" in attribute_list:
            evidence = attribute_list[attribute_list.index("evidence")+1]

        read_list = [getattr(row, "_%s" %x) for x in range(prefix_columns,len(row))]

        ID = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
        reparation_dict[ID] = (proba, evidence, read_list)

    return reparation_dict

def generate_deepribo_dict(deepribo_path):
    """
    create a dictionary containing all important deepribo input
    """

    deepribo_df = pd.read_csv(deepribo_path, sep="\t", comment="#")
    prefix_columns = 9

    deepribo_dict = {}
    for row in deepribo_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        pred_rank = getattr(row, "_5")
        strand = getattr(row, "_6")
        attribute_list = [x for x in re.split('[;=]', getattr(row, "_8")) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")

        pred_value = ""
        if "pred_value" in attribute_list:
            pred_value = attribute_list[attribute_list.index("pred_value")+1]

        if pred_value == "" or float(pred_value) < 0:
            continue

        evidence = ""
        if "evidence" in attribute_list:
            evidence = attribute_list[attribute_list.index("evidence")+1]

        read_list = [getattr(row, "_%s" %x) for x in range(prefix_columns,len(row))]

        ID = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
        deepribo_dict[ID] = (pred_rank, pred_value, evidence, read_list)

    return deepribo_dict

def generate_annotation_dict(annotation_path):
    """
    create dictionary from annotation.
    key : (gene_id, locus_tag, name, gene_name)
    """

    annotation_df = pd.read_csv(annotation_path, sep="\t", comment="#", header=None)
    annotation_dict = {}
    gene_dict = {}
    cds_dict = {}

    for row in annotation_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")
        read_list = [getattr(row, "_%s" %x) for x in range(9,len(row))]

        if ";" in attributes and "=" in attributes:
            attribute_list = [x for x in re.split('[;=]', attributes) if x != ""]
        else:
            attribute_list = [x.replace(";", "") for x in list(csv.reader([attributes], delimiter=' ', quotechar='"'))[0]]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")

        if feature.lower() == "cds":
            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            name = ""
            if "name" in attribute_list:
                name = attribute_list[attribute_list.index("name")+1]

            gene_id = ""
            if "gene_id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("gene_id")+1]
            elif "id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("id")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            cds_dict[new_key] = (gene_id, locus_tag, name, read_list)
        elif feature.lower() in ["gene","pseudogene"]:
            gene_name = ""
            if "name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("name")+1]

            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            gene_dict[new_key] = (gene_name, locus_tag)

    for key in cds_dict.keys():
        gene_name = ""
        gene_id, locus_tag, name, read_list = cds_dict[key]
        if key in gene_dict:
            gene_locus_tag, gene_name = gene_dict[key]
            if locus_tag == "":
                locus_tag = gene_locus_tag
        annotation_dict[key] = (gene_id, locus_tag, name, read_list, gene_name)

    return annotation_dict

def create_excel_file(args):
    # read annotation from file
    annotation_dict = generate_annotation_dict(args.annotation_path)

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

    wildcards = get_unique(wildcards)

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

    conditions = get_unique(conditions)

    xtail_dict, riborex_dict, deepribo_dict, reparation_dict = {}, {}, {}, {}
    if args.xtail_path != "":
        xtail_dict = generate_xtail_dict(args.xtail_path)
    if args.riborex_path != "":
        riborex_dict = generate_riborex_dict(args.riborex_path)
    if args.reads_deepribo != "":
        deepribo_dict = generate_deepribo_dict(args.reads_deepribo)
    if args.reads_reparation != "":
        reparation_dict = generate_reparation_dict(args.reads_reparation)

    #keys_union = list(set(deepribo_dict.keys()) | set(reparation_dict.keys()))

    keys_union = list(set().union(deepribo_dict.keys(), reparation_dict.keys(), annotation_dict.keys()))

    # read gff file
    all_sheet = []

    header = ["Genome", "Start", "Stop", "Strand", "Locus_tag", "Name", "Gene_name", "Length", "Codon_count", "Start_codon", "Stop_codon", "Nucleotide_seq", "Aminoacid_seq"] + [cond + "_TE" for cond in TE_header] + [card + "_rpkm" for card in wildcards] +\
             ["Evidence", "Reparation_probability", "Deepribo_rank", "Deepribo_score", "riborex_pvalue", "riborex_pvalue_adjusted","riborex_log2FC", "xtail_pvalue", "xtail_pvalue_adjusted", "xtail_log2FC"]

    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    for key in keys_union:
        chromosome, mid, strand = key.split(":")
        start, stop = mid.split("-")

        locus_tag, name, gene_name = "","",""
        reparation_probability, deepribo_rank, deepribo_score = 0, 0, 0
        riborex_pvalue, riborex_pvalue_adjusted, riborex_log2FC = 0, 0, 0
        xtail_pvalue, xtail_pvalue_adjusted, xtail_log2FC = 0, 0, 0
        contrast_list = []
        evidence = []

        read_list = []

        gene_id = ""
        if key in annotation_dict:
            gene_id, locus_tag, name, read_list, gene_name = annotation_dict[key]

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
            evidence.extend(evidence_list)

        if key in deepribo_dict:
            deepribo_rank, deepribo_score, deepribo_evidence, read_list = deepribo_dict[key]
            evidence_list = []
            for e in deepribo_evidence.split(" "):
                if not "deepribo" in e:
                    evidence_list.append("deepribo-"+e)
                else:
                    evidence_list.append(e)
            evidence.extend(evidence_list)

        if key in xtail_dict:
            xtail_log2FC, xtail_pvalue, xtail_pvalue_adjusted, contrasts = xtail_dict[key]
            contrast_list.extend(contrasts.split(" "))
        elif gene_id != "" and gene_id in xtail_dict:
            xtail_log2FC, xtail_pvalue, xtail_pvalue_adjusted, contrasts = xtail_dict[gene_id]
            contrast_list.extend(contrasts.split(" "))

        if key in riborex_dict:
            riborex_log2FC, riborex_pvalue, riborex_pvalue_adjusted, contrasts = riborex_dict[key]
            contrast_list.extend(contrasts.split(" "))
        elif gene_id != "" and gene_id in riborex_dict:
            riborex_log2FC, riborex_pvalue, riborex_pvalue_adjusted, contrasts = riborex_dict[gene_id]
            contrast_list.extend(contrasts.split(" "))

        start_codon, stop_codon, nucleotide_seq, aa_seq = get_genome_information(genome_dict[chromosome], int(start)-1, int(stop)-1, strand)

        contrast_list.sort()
        evidence.sort()

        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(calculate_rpkm(total_mapped_dict[(wildcards[idx], chromosome)], val, length))

        TE_list = calculate_TE(rpkm_list, wildcards, conditions)

        evidence = " ".join(evidence)
        result = [chromosome, start, stop, strand, locus_tag, name, gene_name, length, codon_count, start_codon, stop_codon, nucleotide_seq, aa_seq] + TE_list + rpkm_list +\
                 [evidence, reparation_probability, deepribo_rank, deepribo_score, riborex_pvalue, riborex_pvalue_adjusted, riborex_log2FC, xtail_pvalue, xtail_pvalue_adjusted, xtail_log2FC]

        all_sheet.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(all_sheet, columns=[header[x] for x in range(len(header))])

    all_df = all_df.astype({"Start" : "int32", "Stop" : "int32"})
    all_df = all_df.sort_values(by=["Genome", "Start", "Stop"])

    all_df.to_csv(args.output_path.replace(".xlsx", ".tsv"), sep="\t", index=False, quoting=csv.QUOTE_NONE)

    dataframe_dict = { "all" : all_df }

    excel_writer(args, dataframe_dict, wildcards)


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
