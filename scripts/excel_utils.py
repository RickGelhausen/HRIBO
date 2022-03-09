#!/usr/bin/env python

import re
import sys
import pandas as pd
from collections import Counter, OrderedDict


from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

class OrderedCounter(Counter, OrderedDict):
    pass


def get_te_header(wildcards):
    """
    generate the correct TE_header based on the available data
    """
    te_header = []
    te_header_dict = OrderedDict()
    for card in wildcards:
        method, condition, replicate = card.split("-")
        if method == "RIBO":
            if "%s-%s-%s" %("RNA", condition, replicate) in wildcards:
                if ("RIBO", condition) in  te_header_dict:
                    te_header_dict[("RIBO", condition)].append(replicate)
                else:
                    te_header_dict[("RIBO", condition)] = [replicate]
        elif method == "TIS":
            if "%s-%s-%s" %("RNATIS", condition, replicate) in wildcards:
                if ("TIS", condition) in  te_header_dict:
                    te_header_dict[("TIS", condition)].append(replicate)
                else:
                    te_header_dict[("TIS", condition)] = [replicate]
        elif method == "TTS":
            if "%s-%s-%s" %("RNATTS", condition, replicate) in wildcards:
                if ("TTS", condition) in  te_header_dict:
                    te_header_dict[("TTS", condition)].append(replicate)
                else:
                    te_header_dict[("TTS", condition)] = [replicate]

    for key, val in te_header_dict.items():
        method, condition = key
        if len(val) > 1:
            t_header = ["%s-%s-%s" % (method, condition, x) for x in val] + ["%s-%s-avg" % (method, condition)]
        else:
            t_header = ["%s-%s-%s" % (method, condition, x) for x in val]
        te_header.extend(t_header)

    return te_header

def calculate_rpkm(total_mapped, read_count, read_length):
    """
    calculate the rpkm
    """
    if read_length == 0:
        print("read_length: 0 detected! Setting RPKM to 0!")
        return 0
    elif total_mapped == 0:
        print("total_mapped: 0 detected! Setting RPKM to 0!")
        return 0

    return float("%.2f" % ((read_count * 1000000000) / (total_mapped * read_length)))

def get_unique(in_list):
    seen = set()
    seen_add = seen.add
    return sorted([x for x in in_list if not (x in seen or seen_add(x))])

def retrieve_column_information(attributes):
    """
    check for gff2/gff3 format and generate a list of information for the final tables
    [pred_value, name, product, note, evidence, locus_tag, old_locus_tag]
    """

    attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

    if "ORF_type=;" in attributes:
        attribute_list.remove("ORF_type")

    if len(attribute_list) % 2 == 0:
        for i in range(len(attribute_list)):
            if i % 2 == 0:
                attribute_list[i] = attribute_list[i].lower()
    else:
        print(attributes)
        sys.exit("Attributes section of gtf/gff is wrongly formatted!")

    pred_value = ""
    if "pred_value" in attribute_list:
        pred_value = attribute_list[attribute_list.index("pred_value")+1]
    elif "prob" in attribute_list:
        pred_value = attribute_list[attribute_list.index("prob")+1]

    locus_tag = ""
    if "locus_tag" in attribute_list:
        locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

    old_locus_tag = ""
    if "old_locus_tag" in attribute_list:
        old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

    name = ""
    if "name" in attribute_list:
        name = attribute_list[attribute_list.index("name")+1]

    product = ""
    if "product" in attribute_list:
        product = attribute_list[attribute_list.index("product")+1]

    note = ""
    if "note" in attribute_list:
        note = attribute_list[attribute_list.index("note")+1]

    evidence = ""
    if "evidence" in attribute_list:
        evidence = attribute_list[attribute_list.index("evidence")+1]

    return [pred_value, name, product, note, evidence, locus_tag, old_locus_tag]

def get_genome_information(genome, start, stop, strand):
    """
    retrieve the nucleotide sequence and amino acid sequence
    and the start and stop codons
    """
    if strand == "+":
        nucleotide_seq = genome[0][start:stop+1]
        nt_window = genome[0][start-15:start]
    else:
        nucleotide_seq = genome[1][start:stop+1][::-1]
        nt_window = str(Seq(genome[0][stop+1:stop+16]).reverse_complement())

    start_codon = nucleotide_seq[0:3]
    stop_codon = nucleotide_seq[-3:]

    coding_dna = Seq(nucleotide_seq, generic_dna)
    if len(coding_dna) % 3 != 0:
        aa_seq = ""
    else:
        aa_seq = str(coding_dna.translate(table=11,to_stop=True))

    return start_codon, stop_codon, nucleotide_seq, aa_seq, nt_window

def excel_writer(output_path, data_frames, wildcards):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Note", "Aminoacid_seq", "Nucleotide_seq", "Start_codon", "Stop_codon", "Strand", "Codon_count"] + [card + "_rpkm" for card in wildcards]
    writer = pd.ExcelWriter(output_path, engine='xlsxwriter')
    for sheetname, df in data_frames.items():
        df.to_excel(writer, sheet_name=sheetname, index=False)
        worksheet = writer.sheets[sheetname]
        worksheet.freeze_panes(1, 0)
        for idx, col in enumerate(df):
            series = df[col]
            if col in header_only:
                max_len = len(str(series.name)) + 2
            else:
                max_len = max(( series.astype(str).str.len().max(), len(str(series.name)) )) + 1
            #print("Sheet: %s | col: %s | max_len: %s" % (sheetname, col, max_len))
            worksheet.set_column(idx, idx, max_len)
    writer.save()

def te(ribo_count, rna_count):
    """
    calculate the translational efficiency for one entry
    """

    if ribo_count == 0 and rna_count == 0:
        return "NaN"
    elif rna_count == 0:
        return "NaN"
    else:
        return ribo_count / rna_count

def get_avg(t_eff):
    """
    get the final TE list
    """

    valid_count = 0
    sum = 0
    for t in t_eff:
        if t != "NaN":
            valid_count += 1
            sum += t

    if valid_count == 0:
        t_eff.extend(["NaN"])

    else:
        t_eff.extend([sum / valid_count])

    return t_eff

def calculate_te(read_list, wildcards, conditions):
    """
    calculate the translational efficiency
    """
    read_dict = OrderedDict()
    te_dict = OrderedDict()
    for idx in range(len(wildcards)):
        method, condition, replicate = wildcards[idx].split("-")
        key = (method, condition, replicate)
        if key not in read_dict:
            read_dict[key] = read_list[idx]
        else:
            print("warning: multiple equal keys")

    te_list = []
    for key, val in read_dict.items():
        method, condition, replicate = key
        if method == "RIBO":
            if ("RNA", condition, replicate) in read_dict:
                rpkm_ribo = read_dict[key]
                rpkm_rna = read_dict[("RNA", condition, replicate)]
                cur_te = te(rpkm_ribo, rpkm_rna)
                if ("RIBO", condition) in te_dict:
                    te_dict[("RIBO", condition)].append(cur_te)
                else:
                    te_dict[("RIBO", condition)] = [cur_te]

        elif method == "TIS":
            if ("RNATIS", condition, replicate) in read_dict:
                rpkm_ribo = read_dict[key]
                rpkm_rna = read_dict[("RNATIS", condition, replicate)]
                cur_te = te(rpkm_ribo, rpkm_rna)
                if ("TIS", condition) in te_dict:
                    te_dict[("TIS", condition)].append(cur_te)
                else:
                    te_dict[("TIS", condition)] = [cur_te]

        elif method == "TTS":
            if ("RNATTS", condition, replicate) in read_dict:
                rpkm_ribo = read_dict[key]
                rpkm_rna = read_dict[("RNATTS", condition, replicate)]
                cur_te = te(rpkm_ribo, rpkm_rna)
                if ("TTS", condition) in te_dict:
                    te_dict[("TTS", condition)].append(cur_te)
                else:
                    te_dict[("TTS", condition)] = [cur_te]

    te_list = []
    for key, val in te_dict.items():
        if len(val) > 1:
            t_eff = get_avg(val)
        else:
            t_eff = val
        te_list.extend(t_eff)

    return te_list

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
        contrast = getattr(row, "contrast").split("_")[1]

        riborex_dict[(gene_id, contrast)] = (log2fc, pvalue, pvalue_adj)

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
        contrast = getattr(row, "contrast").split("_")[1]

        xtail_dict[(gene_id, contrast)] = (log2fc, pvalue, pvalue_adj)

    return xtail_dict


def generate_reparation_dict(reparation_path):
    """
    create a dictionary containing all important reparation input
    """

    reparation_df = pd.read_csv(reparation_path, header=None, sep="\t", comment="#")
    prefix_columns = 9

    reparation_dict = {}
    for row in reparation_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attribute_list = [x.strip(" ") for x in re.split('[;=]', getattr(row, "_8")) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attribute_list)
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

    deepribo_df = pd.read_csv(deepribo_path, header=None, sep="\t", comment="#")
    prefix_columns = 9

    deepribo_dict = {}
    for row in deepribo_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        pred_rank = getattr(row, "_5")
        strand = getattr(row, "_6")
        attribute_list = [x.strip(" ") for x in re.split('[;=]', getattr(row, "_8")) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attribute_list)
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")

        pred_value = ""
        if "pred_value" in attribute_list:
            pred_value = attribute_list[attribute_list.index("pred_value")+1]

        if pred_value == "":# or float(pred_value) < 0:
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
    annotation_meta_dict = {}

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

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attribute_list)
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")


        if feature.lower() == "cds":
            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            name = ""
            if "name" in attribute_list:
                name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                name = attribute_list[attribute_list.index("gene_name")+1]

            gene_id = ""
            if "gene_id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("gene_id")+1]
            elif "id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("id")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            cds_dict[new_key] = (gene_id, locus_tag, name, read_list, old_locus_tag)
        elif feature.lower() in ["gene","pseudogene"]:
            gene_name = ""
            if "name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("gene_name")+1]

            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]
            elif "gene_id" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("gene_id")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            gene_dict[new_key] = (gene_name, locus_tag, old_locus_tag)


    for key in cds_dict.keys():
        gene_name = ""
        gene_id, locus_tag, name, read_list, old_locus_tag = cds_dict[key]

        if key in gene_dict:
            gene_name, gene_locus_tag, gene_old_locus_tag = gene_dict[key]

            if locus_tag == "":
                locus_tag = gene_locus_tag
            if old_locus_tag == "":
                old_locus_tag = gene_old_locus_tag

        annotation_meta_dict[key] = (gene_id, locus_tag, name, gene_name, old_locus_tag, read_list)

    return annotation_meta_dict

def generate_non_cds_dict(annotation_path):
    """
    create dictionary from annotation ignoring cds.
    key : (gene_id, locus_tag, name, gene_name)
    """

    gene_dict = {}
    non_cds_dict = {}

    annotation_meta_dict = {}
    annotation_df = pd.read_csv(annotation_path, sep="\t", comment="#", header=None)

    for row in annotation_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")
        read_list = [getattr(row, "_%s" %x) for x in range(9,len(row))]

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attribute_list)
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")


        if feature.lower() not in ["cds", "gene", "pseudogene", "exon"]:
            locus_tag = None
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            old_locus_tag = None
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            name = None
            if "name" in attribute_list:
                name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                name = attribute_list[attribute_list.index("gene_name")+1]

            gene_id = None
            if "gene_id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("gene_id")+1]
            elif "id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("id")+1]

            key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            non_cds_dict[key] = (feature, gene_id, locus_tag, old_locus_tag, name, read_list)

        elif feature.lower() in ["gene","pseudogene"]:
            gene_name = ""
            if "name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("gene_name")+1]

            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]
            elif "gene_id" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("gene_id")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            gene_dict[new_key] = (gene_name, locus_tag, old_locus_tag)


    for key in non_cds_dict.keys():
        gene_name = ""
        feature, gene_id, locus_tag, old_locus_tag, name, read_list = non_cds_dict[key]

        if key in gene_dict:
            gene_name, gene_locus_tag, gene_old_locus_tag = gene_dict[key]

            if locus_tag == "":
                locus_tag = gene_locus_tag
            if old_locus_tag == "":
                old_locus_tag = gene_old_locus_tag

        annotation_meta_dict[key] = (feature, gene_id, locus_tag, name, gene_name, old_locus_tag, read_list)

    return annotation_meta_dict
