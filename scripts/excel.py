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
    return "%.2f" % ((read_count * 1000000000) / (total_mapped * read_length))

def get_genome_information(genome, start, stop, strand):
    # retrieve the nucleotide sequence with correct strand
    if strand == "+":
        nucleotide_seq = genome[0][start:stop+1]
    else:
        nucleotide_seq = genome[1][start:stop+1]

    # get start and stop codons
    start_codon = nucleotide_seq[0:3]
    stop_codon = nucleotide_seq[-3:]

    # translate nucleotide_seq to aminoacid sequence
    coding_dna = Seq(nucleotide_seq, generic_dna)
    aa_seq = str(coding_dna.translate(table=11))
    return start_codon, stop_codon, nucleotide_seq, aa_seq

def parse_orfs(args):
    # read the genome file
    genome_file = SeqIO.parse(args.genome, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq),str(entry.seq)[::-1])

    # get the total mapped reads for each bam file
    total_mapped_dict = {}
    with open(args.total_mapped, "r") as f:
        total = f.readlines()

    for line in total:
        key, value = line.strip().split("\t")
        total_mapped_dict[key] = int(value)

    # read the comment containing the wildcards
    with open(args.reads, "r") as f:
        wildcards = f.readline()

    wildcards = wildcards.replace("# ", "").split()
    wildcards = [os.path.splitext(os.path.basename(card))[0] for card in wildcards]

    read_df = pd.read_csv(args.reads, comment="#", header=None, sep="\t")

    # id + start + stop + strand + length + rest
    column_count = len(read_df.columns) + 8                      # ADD HERE

    name_list = ["s%s" % str(x) for x in range(column_count)]
    nTuple = collections.namedtuple('Pandas', name_list)
    # read gff file
    rows = []
    header = ["orfID", "start", "stop", "strand", "length"] + [card + "_rpkm" for card in wildcards] + ["evidence", "annotated", "name", "ORF_type", "start_codon", "stop_codon", "nucleotide_seq", "aminoacid_seq"] # ADD HERE
    rows.append(nTuple(*header))
    for row in read_df.itertuples(index=False, name='Pandas'):
        accession = getattr(row, "_0")
        start = getattr(row, "_1")
        stop = getattr(row, "_2")
        strand = getattr(row, "_3")
        attributes = getattr(row, "_4")

        start_codon, stop_codon, nucleotide_seq, aa_seq = get_genome_information(genome_dict[accession], start, stop, strand)

        attribute_list = re.split('[;=]', attributes)
        id = attribute_list[attribute_list.index("ID")+1]

        annotated="NA"
        if "annotated" in attribute_list:
            annotated = attribute_list[attribute_list.index("annotated")+1]

        name="NA"
        if "name" in attribute_list:
            name = attribute_list[attribute_list.index("name")+1]

        evidence = attribute_list[attribute_list.index("Evidence")+1]
        orftype = attribute_list[attribute_list.index("ORF_type")+1]
        if orftype == "":
            orftype = "NA"

        length = stop - start + 1

        read_list = [getattr(row, "_%s" %x) for x in range(5,len(row))]

        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(calculate_rpkm(total_mapped_dict[wildcards[idx]], val, length))

        result = [id, start, stop, strand, length] + rpkm_list + [evidence, annotated, name, orftype, start_codon, stop_codon, nucleotide_seq, aa_seq] # ADD HERE
        rows.append(nTuple(*result))

    excel_df = pd.DataFrame.from_records(rows, columns=[x for x in range(column_count)])

    excel_df.to_excel(args.output, sheet_name=args.sheet_name)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel files containing: \
                                                id, start, stop, orflength, potential RBS, rpkm')
    #parser.add_argument("-a", "--annotation_combine", action="store", dest="combined", required=True, help= "orf gff file.")
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True, help= "reference genome")
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped", required=True\
                                                    , help= "file containing the total mapped reads for all alignment files.")
    parser.add_argument("-r", "--mapped_reads", action="store", dest="reads", required=True, help= "file containing the individual read counts")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    parser.add_argument("-s", "--sheet_name", action="store", dest="sheet_name", help= "name of the excel sheet", default="Sheet 1")

    args = parser.parse_args()

    parse_orfs(args)

if __name__ == '__main__':
    main()
