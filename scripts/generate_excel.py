#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import collections

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

def calculate_rpkm(total_mapped, read_count, read_length):
    """
    calculate the rpkm
    """
    return float("%.2f" % ((read_count * 1000000000) / (total_mapped * read_length)))

def retrieve_column_information(attributes):
    """
    check for gff2/gff3 format and generate a list of information for the final tables
    [locus_tag, name, product, note]
    """

    if ";" in attributes and "=" in attributes:
        attribute_list =  [x for x in re.split('[;=]', attributes)]
    else:
        attribute_list = [x.replace("\"", "") for x in re.split('[; ]', attributes) if x != ""]

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

    product = ""
    if "product" in attribute_list:
        product = attribute_list[attribute_list.index("product")+1]

    note = ""
    if "note" in attribute_list:
        note = attribute_list[attribute_list.index("note")+1]

    evidence = ""
    if "evidence" in attribute_list:
        evidence = attribute_list[attribute_list.index("evidence")+1]

    orftype = ""
    if "orf_type" in attribute_list:
        orftype = attribute_list[attribute_list.index("orf_type")+1]

    return [locus_tag, name, product, note, evidence, orftype]

def get_genome_information(genome, start, stop, strand):
    # retrieve the nucleotide sequence with correct strand
    if strand == "+":
        nucleotide_seq = genome[0][start:stop+1]
    else:
        nucleotide_seq = genome[1][start:stop+1][::-1]

    # get start and stop codons
    start_codon = nucleotide_seq[0:3]
    stop_codon = nucleotide_seq[-3:]

    # translate nucleotide_seq to aminoacid sequence
    coding_dna = Seq(nucleotide_seq, generic_dna)
    aa_seq = str(coding_dna.translate(table=11,to_stop=True))
    return start_codon, stop_codon, nucleotide_seq, aa_seq

def excel_writer(args, data_frames, wildcards):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Note", "Aminoacid_seq", "Nucleotide_seq", "Start_codon", "Stop_codon", "Strand", "Codon_count"] + [card + "_rpkm" for card in wildcards]
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

def parse_orfs(args):
    # read the genome file
    genome_file = SeqIO.parse(args.genome, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))
        with open("test-rev.txt", "w") as f:
            f.write(str(entry.seq.complement()))

    # get the total mapped reads for each bam file
    total_mapped_dict = {}
    with open(args.total_mapped, "r") as f:
        total = f.readlines()

    for line in total:
        wildcard, reference_name, value = line.strip().split("\t")
        total_mapped_dict[(wildcard, reference_name)] = int(value)

    # read the comment containing the wildcards
    with open(args.reads, "r") as f:
        wildcards = f.readline()

    wildcards = wildcards.replace("# ", "").split()
    wildcards = [os.path.splitext(os.path.basename(card))[0] for card in wildcards]

    #read bed file
    read_df = pd.read_csv(args.reads, comment="#", header=None, sep="\t")

    # read gff file
    main_sheet = []
    cds_sheet = []
    gene_sheet = []
    region_sheet = []
    rRNA_sheet = []
    tRNA_sheet = []
    misc_sheet = []

    header = ["Genome", "Source", "Feature", "Start", "Stop", "Strand", "Locus_tag", "Name", "Length", "Codon_count"] + [card + "_rpkm" for card in wildcards] + ["Evidence", "ORF_type", "Start_codon", "Stop_codon", "Nucleotide_seq", "Aminoacid_seq",  "Product", "Note"]
    prefix_columns = len(read_df.columns) - len(wildcards)
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    for row in read_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        start = getattr(row, "_1")
        stop = getattr(row, "_2")
        strand = getattr(row, "_5")

        attributes = getattr(row, "_3")
        source = getattr(row, "_6")
        feature = getattr(row, "_7")

        start_codon, stop_codon, nucleotide_seq, aa_seq = get_genome_information(genome_dict[reference_name], start-1, stop-1, strand)
        column_info = retrieve_column_information(attributes)

        length = stop - start + 1
        codon_count = int(length / 3)

        read_list = [getattr(row, "_%s" %x) for x in range(prefix_columns,len(row))]
        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(calculate_rpkm(total_mapped_dict[(wildcards[idx], reference_name)], val, length))

        result = [reference_name, source, feature, start, stop, strand, column_info[0], column_info[1], length, codon_count] + rpkm_list + [column_info[4], column_info[5], start_codon, stop_codon, nucleotide_seq, aa_seq, column_info[2], column_info[3]]
        if feature.lower() != "region":
            main_sheet.append(nTuple(*result))

        if feature.lower() == "cds":
            cds_sheet.append(nTuple(*result))
        elif feature.lower() == "gene":
            gene_sheet.append(nTuple(*result))
        elif feature.lower() == "region":
            region_sheet.append(nTuple(*result))
        elif feature.lower() == "rrna":
            rRNA_sheet.append(nTuple(*result))
        elif feature.lower() == "trna":
            tRNA_sheet.append(nTuple(*result))
        else:
            misc_sheet.append(nTuple(*result))

    main_df = pd.DataFrame.from_records(main_sheet, columns=[header[x] for x in range(len(header))])
    cds_df = pd.DataFrame.from_records(cds_sheet, columns=[header[x] for x in range(len(header))])
    gene_df = pd.DataFrame.from_records(gene_sheet, columns=[header[x] for x in range(len(header))])
    region_df = pd.DataFrame.from_records(region_sheet, columns=[header[x] for x in range(len(header))])
    rRNA_df = pd.DataFrame.from_records(rRNA_sheet, columns=[header[x] for x in range(len(header))])
    tRNA_df = pd.DataFrame.from_records(tRNA_sheet, columns=[header[x] for x in range(len(header))])
    misc_df = pd.DataFrame.from_records(misc_sheet, columns=[header[x] for x in range(len(header))])

    dataframe_dict = { "Main" : main_df, "CDS" : cds_df, "gene" : gene_df, "region" : region_df, "rRNA" : rRNA_df, "tRNA" : tRNA_df, "miscellaneous" : misc_df }

    excel_writer(args, dataframe_dict, wildcards)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel files containing: \
                                                id, start, stop, orflength, potential RBS, rpkm')
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True, help= "reference genome")
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped", required=True\
                                                    , help= "file containing the total mapped reads for all alignment files.")
    parser.add_argument("-r", "--mapped_reads", action="store", dest="reads", required=True, help= "file containing the individual read counts")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    args = parser.parse_args()

    parse_orfs(args)

if __name__ == '__main__':
    main()
