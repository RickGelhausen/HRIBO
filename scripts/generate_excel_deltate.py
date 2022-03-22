#!/usr/bin/env python
import argparse
import sys
import pandas as pd
import collections

import excel_utils as eu

from Bio import SeqIO

def create_combined_dict(ribo_df, rna_df, te_df):
    """
    Merge into one dataframe
    """

    result_dict = {}
    for row in ribo_df.itertuples(index=False):
        identifier = getattr(row, "Identifier")
        if identifier in result_dict:
            _, tmp_rna, tmp_te = result_dict[identifier]
            result_dict[identifier] = [row, tmp_rna, tmp_te]
        else:
            result_dict[identifier] = [row, None, None]

    for row in rna_df.itertuples(index=False):
        identifier = getattr(row, "Identifier")
        if identifier in result_dict:
            tmp_ribo, _, tmp_te = result_dict[identifier]
            result_dict[identifier] = [tmp_ribo, row, tmp_te]
        else:
            result_dict[identifier] = [None, row, None]

    for row in te_df.itertuples(index=False):
        identifier = getattr(row, "Identifier")
        if identifier in result_dict:
            tmp_ribo, tmp_rna, _ = result_dict[identifier]
            result_dict[identifier] = [tmp_ribo, tmp_rna, row]
        else:
            result_dict[identifier] = [None, None, row]

    return result_dict

def deltate_output(args):
    # read the genome file
    genome_file = SeqIO.parse(args.genome, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    annotation_dict = eu.annotation_to_dict(args.annotation_file)

    ribo_expr_df = pd.read_csv(args.input_ribo, sep="\t", comment="#")
    ribo_expr_df.index.name = "Identifier"
    ribo_expr_df = ribo_expr_df.reset_index(level=["Identifier"])

    rna_expr_df = pd.read_csv(args.input_rna, sep="\t", comment="#")
    rna_expr_df.index.name = "Identifier"
    rna_expr_df = rna_expr_df.reset_index(level=["Identifier"])

    te_expr_df = pd.read_csv(args.input_te, sep="\t", comment="#")
    te_expr_df.index.name = "Identifier"
    te_expr_df = te_expr_df.reset_index(level=["Identifier"])

    combined_dict = create_combined_dict(ribo_expr_df, rna_expr_df, te_expr_df)

    all_sheet = []
    header = ["Genome", "Start", "Stop", "Strand", "Locus_tag", "Old_locus_tag", "Identifier", "Name", \
              "RIBO_baseMean", "RIBO_log2FoldChange", "RIBO_lfcSE", "RIBO_pvalue", "RIBO_padj", \
              "RNA_baseMean", "RNA_log2FoldChange", "RNA_lfcSE", "RNA_pvalue", "RNA_padj", \
              "TE_baseMean", "TE_log2FoldChange", "TE_lfcSE", "TE_stat", "TE_pvalue", "TE_padj", \
              "Length", "Codon_count", "Start_codon", "Stop_codon", "Nucleotide_seq", "Aminoacid_seq"]
    name_list = [f"s{x}" for x in range(len(header))]

    nTuple = collections.namedtuple('Pandas', name_list)

    for unique_id, rows in combined_dict.items():
        if unique_id in annotation_dict:
            chromosome, start, stop, strand, gene_name, locus_tag, old_locus_tag = annotation_dict[unique_id]
        else:
            if ":" in unique_id and "-" in unique_id:
                chromosome, sec, strand = unique_id.split(":")
                start, stop = sec.split("-")
                gene_name = ""
                locus_tag, old_locus_tag = "", ""
            else:
                sys.exit("Error... ID is not novel and not in the annotation!")

        if None in rows:
            continue

        ribo_list = list(rows[0])[1:]
        rna_list = list(rows[1])[1:]
        te_list = list(rows[2])[1:]

        start = int(start)
        stop = int(stop)
        length = stop - start + 1
        codon_count = int(length / 3)
        start_codon, stop_codon, nucleotide_seq, aa_seq, _ = eu.get_genome_information(genome_dict[chromosome], start-1, stop-1, strand)

        result = [chromosome, start, stop, strand, locus_tag, old_locus_tag, unique_id, gene_name] \
               + ribo_list + rna_list + te_list \
               + [length, codon_count, start_codon, stop_codon, nucleotide_seq, aa_seq]

        all_sheet.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(all_sheet, columns=[header[x] for x in range(len(header))])
    all_df = all_df.sort_values(by=["TE_padj", "Genome", "Start", "Stop", "Strand"])
    significant_df = all_df[all_df["TE_padj"] <= 0.05]
    dataframe_dict = {"all" : all_df, "significant" : significant_df}

    eu.excel_writer(args.output, dataframe_dict, [])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel files from riborex output')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_file", required=True, help= "annotation file")
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True, help= "reference genome")
    parser.add_argument("-i", "--delta_ribo", action="store", dest="input_ribo", required=True, help= "input txt file for ribo")
    parser.add_argument("-r", "--delta_rna", action="store", dest="input_rna", required=True, help= "input txt file for rna")
    parser.add_argument("-t", "--delta_te", action="store", dest="input_te", required=True, help= "input txt file for te")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    args = parser.parse_args()

    deltate_output(args)

if __name__ == '__main__':
    main()
