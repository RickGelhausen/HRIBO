#!/usr/bin/env python
import argparse
import sys
import pandas as pd
import collections

import excel_utils as eu

from Bio import SeqIO

def riborex_output(args):
    # read the genome file
    genome_file = SeqIO.parse(args.genome, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    annotation_dict = eu.generate_annotation_dict(args.annotation_file)

    diff_expr_df = pd.read_csv(args.input_csv, sep=",", comment="#")

    all_sheet = []
    header = ["Genome", "Start", "Stop", "Strand", "Locus_tag", "Old_locus_tag", "ID", "Name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Length", "Codon_count", "Start_codon", "Stop_codon", "Nucleotide_seq", "Aminoacid_seq"]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    for row in diff_expr_df.itertuples(index=False, name='Pandas'):
        unique_id = getattr(row, "_0")

        if unique_id in annotation_dict:
            chromosome, start, stop, strand, attributes, locus_tag, old_locus_tag = annotation_dict[unique_id]
            _, gene_name, _, _, _, locus_tag, old_locus_tag = eu.retrieve_column_information(attributes)

        else:
            if ":" in unique_id and "-" in unique_id:
                chromosome, sec, strand = unique_id.split(":")
                start, stop = sec.split("-")
                gene_name = ""
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
        start_codon, stop_codon, nucleotide_seq, aa_seq, _ = eu.get_genome_information(genome_dict[chromosome], start-1, stop-1, strand)

        result = [chromosome, start, stop, strand, locus_tag, old_locus_tag, unique_id, gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, length, codon_count, start_codon, stop_codon, nucleotide_seq, aa_seq]

        all_sheet.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(all_sheet, columns=[header[x] for x in range(len(header))])
    all_df = all_df.sort_values(by=["padj", "Genome", "Start", "Stop", "Strand"])
    significant_df = all_df[all_df["padj"] <= 0.05]
    dataframe_dict = {"all" : all_df, "significant" : significant_df}

    eu.excel_writer(args.output, dataframe_dict, [])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel files from riborex output')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_file", required=True, help= "annotation file")
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True, help= "reference genome")
    parser.add_argument("-i", "--input", action="store", dest="input_csv", required=True, help= "input csv file")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    args = parser.parse_args()

    riborex_output(args)

if __name__ == '__main__':
    main()