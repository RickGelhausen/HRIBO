#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections
import csv

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna


def excel_writer(args, data_frames):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Note", "Aminoacid_seq", "Nucleotide_seq", "Start_codon", "Stop_codon", "Strand", "Codon_count"]
    writer = pd.ExcelWriter(args.output, engine='xlsxwriter')
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



def retrieve_column_information(attributes):
    """
    check for gff2/gff3 format and generate a list of information for the final tables
    [locus_tag, name, product, note]
    """

    attribute_list = [x for x in re.split('[;=]', attributes) if x != ""]

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

    old_locus_tag = ""
    if "old_locus_tag" in attribute_list:
        old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

    name = ""
    if "name" in attribute_list:
        name = attribute_list[attribute_list.index("name")+1]

    return [name]

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
        aa_seq = str(coding_dna.translate(table=11,to_stop=False))
    return start_codon, stop_codon, nucleotide_seq, aa_seq

def annotation_to_dict(annotation_file):
    """
    read the annotation file into a dictionary
    { ID : (chromosome,start,stop,strand,attributes) }
    """

    annotation_df = pd.read_csv(annotation_file, sep="\t", comment="#", header=None)

    parent_dict = {}
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        if getattr(row, "_2").lower() in ["gene","pseudogene"]:
            chromosome = getattr(row, "_0")
            start = getattr(row, "_3")
            stop = getattr(row, "_4")
            strand = getattr(row, "_6")
            attributes = getattr(row, "_8")

            attribute_list = [x for x in re.split('[;=]', attributes) if x != ""]

            id = attribute_list[attribute_list.index("ID") + 1]

            locus_tag = ""
            if "locus_tag" in attributes.lower():
                locus_tag = attribute_list[attribute_list.index("locus_tag") + 1]

            old_locus_tag = ""
            if "old_locus_tag" in attributes.lower():
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag") + 1]

            parent_dict[id] = (locus_tag, old_locus_tag)

    annotation_dict = {}
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        if getattr(row, "_2").lower() == "cds":
            chromosome = getattr(row, "_0")
            start = getattr(row, "_3")
            stop = getattr(row, "_4")
            strand = getattr(row, "_6")
            attributes = getattr(row, "_8")


            attribute_list = [x for x in re.split('[;=]', attributes) if x != ""]

            parent = ""
            if "parent" in attributes.lower():
                parent = attribute_list[attribute_list.index("Parent") + 1]

            id = attribute_list[attribute_list.index("ID") + 1]

            if parent == "" or parent not in parent_dict:
                annotation_dict[id] = (chromosome, start, stop, strand, attributes, "", "")
            else:
                annotation_dict[id] = (chromosome, start, stop, strand, attributes, parent_dict[parent][0], parent_dict[parent][1])


    return annotation_dict

def riborex_output(args):
    # read the genome file
    genome_file = SeqIO.parse(args.genome, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    annotation_dict = annotation_to_dict(args.annotation_file)

    diff_expr_df = pd.read_csv(args.input_csv, sep=",", comment="#")

    all_sheet = []
    header = ["Genome", "Start", "Stop", "Strand", "Locus_tag", "Old_locus_tag", "ID", "Name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Length", "Codon_count", "Start_codon", "Stop_codon", "Nucleotide_seq", "Aminoacid_seq"]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    for row in diff_expr_df.itertuples(index=False, name='Pandas'):
        cds_id = getattr(row, "_0")

        if cds_id in annotation_dict:
            chromosome, start, stop, strand, attributes, locus_tag, old_locus_tag = annotation_dict[cds_id]
            column_info = retrieve_column_information(attributes)

        else:
            if ":" in cds_id and "-" in cds_id:
                chromosome, sec, strand = cds_id.split(":")
                start, stop = sec.split("-")
                column_info = ["","",""]
                locus_tag, old_locus_tag = "", ""
            else:
                sys.exit("Error... ID is not novel and not in the annotation!")

        baseMean = getattr(row, "baseMean")
        log2FoldChange = getattr(row, "log2FoldChange")
        lfcSE = getattr(row, "lfcSE")
        stat = getattr(row, "stat")
        pvalue = getattr(row, "pvalue")
        padj = getattr(row, "padj")

        start = int(start)
        stop = int(stop)
        length = stop - start + 1
        codon_count = int(length / 3)
        start_codon, stop_codon, nucleotide_seq, aa_seq = get_genome_information(genome_dict[chromosome], start-1, stop-1, strand)

        result = [chromosome, start, stop, strand, locus_tag, old_locus_tag, cds_id, column_info[0], baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, length, codon_count, start_codon, stop_codon, nucleotide_seq, aa_seq]

        all_sheet.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(all_sheet, columns=[header[x] for x in range(len(header))])

    dataframe_dict = {"all" : all_df }

    excel_writer(args, dataframe_dict)


def xtail_output(args):
    # read the genome file
    genome_file = SeqIO.parse(args.genome, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    annotation_dict = annotation_to_dict(args.annotation_file)

    diff_expr_df = pd.read_csv(args.input_csv, sep=",", comment="#")

    all_sheet = []
    header = ["Genome", "Start", "Stop", "Strand", "Locus_tag", "Old_locus_tag", "ID", "Name", "mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1", "log2FC_TE_v2", "pvalue_v2", "log2FC_TE_final", "pvalue_final", "pvalue.adjust", "Length", "Codon_count", "Start_codon", "Stop_codon", "Nucleotide_seq", "Aminoacid_seq"]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    for row in diff_expr_df.itertuples(index=False, name='Pandas'):
        cds_id = getattr(row, "_0")

        if cds_id in annotation_dict:
            chromosome, start, stop, strand, attributes, locus_tag, old_locus_tag = annotation_dict[cds_id]
            column_info = retrieve_column_information(attributes)

        else:
            if ":" in cds_id and "-" in cds_id:
                chromosome, sec, strand = cds_id.split(":")
                start, stop = sec.split("-")
                column_info = ["","",""]
                locus_tag, old_locus_tag = "", ""
            else:
                sys.exit("Error... ID is not novel and not in the annotation!")

        mRNA_log2FC = getattr(row, "mRNA_log2FC")
        RPF_log2FC = getattr(row, "RPF_log2FC")
        log2FC_TE_v1 = getattr(row, "log2FC_TE_v1")
        pvalue_v1 = getattr(row, "pvalue_v1")
        log2FC_TE_v2 = getattr(row, "log2FC_TE_v2")
        pvalue_v2 = getattr(row, "pvalue_v2")
        log2FC_TE_final = getattr(row, "log2FC_TE_final")
        pvalue_final = getattr(row, "pvalue_final")
        pvalue_adjust = getattr(row, "_9")

        start = int(start)
        stop = int(stop)
        length = stop - start + 1
        codon_count = int(length / 3)
        start_codon, stop_codon, nucleotide_seq, aa_seq = get_genome_information(genome_dict[chromosome], start-1, stop-1, strand)

        result = [chromosome, start, stop, strand, locus_tag, old_locus_tag, cds_id, column_info[0], mRNA_log2FC, RPF_log2FC,log2FC_TE_v1, pvalue_v1,log2FC_TE_v2,pvalue_v2,log2FC_TE_final,pvalue_final,pvalue_adjust, length, codon_count, start_codon, stop_codon, nucleotide_seq, aa_seq]

        all_sheet.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(all_sheet, columns=[header[x] for x in range(len(header))])

    dataframe_dict = {"all" : all_df }

    excel_writer(args, dataframe_dict)


def create_excel_file(args):
    """
    write excel file depending on the tool in use
    """
    if args.tool == "riborex":
        riborex_output(args)
    else:
        xtail_output(args)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel files from xtail/riborex output')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_file", required=True, help= "annotation file")
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True, help= "reference genome")
    parser.add_argument("-t", "--tool", action="store", dest="tool", required=True, help= "current tool riborex/xtail")
    parser.add_argument("-i", "--input", action="store", dest="input_csv", required=True, help= "input csv file")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    args = parser.parse_args()

    create_excel_file(args)

if __name__ == '__main__':
    main()
