#!/usr/bin/env python
import argparse
import pandas as pd
import collections
from collections import Counter, OrderedDict

from Bio import SeqIO

import excel_utils as eu

class OrderedCounter(Counter, OrderedDict):
    pass


def create_excel_file(args):
    # read the genome file
    genome_file = SeqIO.parse(args.genome, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    # get the total mapped reads for each bam file
    total_mapped_dict = {}
    with open(args.total_mapped, "r") as f:
        total = f.readlines()

    wildcards = []
    for line in total:
        wildcard, chromosome, value = line.strip().split("\t")
        total_mapped_dict[(wildcard, chromosome)] = int(value)
        wildcards.append(wildcard)

    wildcards = eu.get_unique(wildcards)

    te_header = eu.get_te_header(wildcards)

    conditions = []
    for card in wildcards:
        conditions.append(card.split("-")[1])

    conditions = eu.get_unique(conditions)

    #read bed file
    read_df = pd.read_csv(args.reads, comment="#", header=None, sep="\t")

    # read gff file
    cds_sheet = []

    header = ["Identifier", "Genome", "Source", "Feature", "Start", "Stop", "Strand", "Pred_probability", "Locus_tag", "Old_locus_tag", "Name", "Length", "Codon_count"] + [cond + "_TE" for cond in te_header] + [card + "_rpkm" for card in wildcards] + ["Evidence", "Start_codon", "Stop_codon", "15nt upstream", "Nucleotide_seq", "Aminoacid_seq"]
    prefix_columns = len(read_df.columns) - len(wildcards)
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    for row in read_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        source = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")

        start_codon, stop_codon, nucleotide_seq, aa_seq, nt_window = "", "", "", "", ""
        if chromosome in genome_dict:
            start_codon, stop_codon, nucleotide_seq, aa_seq, nt_window = eu.get_genome_information(genome_dict[chromosome], start-1, stop-1, strand)
        pred_value, name, product, note, evidence, locus_tag, old_locus_tag = eu.retrieve_column_information(attributes)

        length = stop - start + 1
        codon_count = int(length / 3)

        read_list = [getattr(row, "_%s" %x) for x in range(prefix_columns,len(row))]
        rpkm_list = []
        for idx, val in enumerate(read_list):
            if (wildcards[idx], chromosome) not in total_mapped_dict:
                rpkm_list.append(0)
            else:
                rpkm_list.append(eu.calculate_rpkm(total_mapped_dict[(wildcards[idx], chromosome)], val, length))

        te_list = eu.calculate_te(rpkm_list, wildcards, conditions)

        identifier = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
        result = [identifier, chromosome, "reparation", feature, start, stop, strand, pred_value, locus_tag, old_locus_tag, name, length, codon_count] + te_list + rpkm_list + [evidence, start_codon, stop_codon, nt_window, nucleotide_seq, aa_seq]


        cds_sheet.append(nTuple(*result))

    cds_df = pd.DataFrame.from_records(cds_sheet, columns=[header[x] for x in range(len(header))])
    cds_df = cds_df.sort_values(by=["Genome", "Start", "Stop"])

    dataframe_dict = { "CDS" : cds_df }

    eu.excel_writer(args.output_path, dataframe_dict, wildcards)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel files containing: \
                                                id, start, stop, orflength, rpkm')
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True, help= "reference genome")
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped", required=True\
                                                    , help= "file containing the total mapped reads for all alignment files.")
    parser.add_argument("-r", "--mapped_reads", action="store", dest="reads", required=True, help= "file containing the individual read counts")
    parser.add_argument("-o", "--xlsx", action="store", dest="output_path", required=True, help= "output xlsx file")
    args = parser.parse_args()

    create_excel_file(args)

if __name__ == '__main__':
    main()
