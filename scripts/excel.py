#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import collections

# from Bio.Seq import Seq
# from Bio import SeqIO

def calculate_rpkm(total_mapped, read_count, read_length):
    """
    calculate the rpkm
    """
    return "%.2f" % ((read_count * 1000000000) / (total_mapped * read_length))

def parse_orfs(args):
    # read the genome file
    # genome_file = SeqIO.parse(args.genome, "fasta")
    # genome_dict = dict()
    # for entry in genome_file:
    #     genome_dict[str(entry.id)] = str(entry.seq)

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
    column_count = len(read_df.columns) + 2

    name_list = ["s%s" % str(x) for x in range(column_count)]
    nTuple = collections.namedtuple('Pandas', name_list)
    # read gff file
    rows = []
    header = ["orfID", "start", "stop", "strand", "length"] + wildcards + ["annotated"]
    rows.append(nTuple(*header))
    for row in read_df.itertuples(index=False, name='Pandas'):
        start = getattr(row, "_1")
        stop = getattr(row, "_2")
        strand = getattr(row, "_3")
        attributes = getattr(row, "_4")

        attribute_list = re.split('[;=]', attributes)
        id = attribute_list[attribute_list.index("ID")+1]

        annotated="NA"
        if "annotated" in attribute_list:
            annotated = attribute_list[attribute_list.index("annotated")+1]

        evidence = attribute_list[attribute_list.index("Evidence")+1]
        length = stop - start + 1

        read_list = [getattr(row, "_%s" %x) for x in range(5,len(row))]

        rpkm_list = []
        for idx, val in enumerate(read_list):
            rpkm_list.append(calculate_rpkm(total_mapped_dict[wildcards[idx]], val, length))

        result = [id, start, stop, strand, length] + rpkm_list + [evidence, annotated]
        rows.append(nTuple(*result))

    excel_df = pd.DataFrame.from_records(rows, columns=[x for x in range(column_count)])

    excel_df.to_excel(args.output, sheet_name=args.sheet_name)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel files containing \
                                                id, start, stop, orflength, potential RBS, rpkm')
    #parser.add_argument("-a", "--annotation_combine", action="store", dest="combined", required=True, help= "orf gff file.")
    #parser.add_argument("-g", "--genome", action="store", dest="genome", required=True, help= "reference genome")
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped", required=True\
                                                    , help= "file containing the total mapped reads for all alignment files.")
    parser.add_argument("-r", "--mapped_reads", action="store", dest="reads", required=True, help= "file containing the individual read counts")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    parser.add_argument("-s", "--sheet_name", action="store", dest="sheet_name", help= "name of the excel sheet", default="Sheet 1")

    args = parser.parse_args()

    parse_orfs(args)

if __name__ == '__main__':
    main()
