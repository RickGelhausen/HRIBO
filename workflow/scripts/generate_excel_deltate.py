#!/usr/bin/env python
"""Turn the deltaTE differential expression results into a spreadsheet.

deltaTE reports three separate DESeq2 tables -- ribosome occupancy, RNA
abundance and their ratio -- which are joined here on the feature identifier.
The genome lookup and identifier resolution are shared with the riborex and
xtail tables in excel_utils.
"""

import argparse

import pandas as pd

import excel_utils as eu

# The RIBO and RNA tables carry no Wald statistic; only the TE table does.
# These are DESeq2's own column names as they appear in the input files.
RIBO_RNA_COLUMNS = ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]
TE_COLUMNS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# DESeq2's spelling -> the workflow-wide header name.
HEADER_NAMES = {
    "log2FoldChange": "log2FC",
    "lfcSE": "log2FC_SE",
    "padj": "pvalue_adjusted",
}


def header_name(prefix, column):
    return f"{prefix}_{HEADER_NAMES.get(column, column)}"

# Each split is (sheet name, fold change column, adjusted p-value column).
SPLITS = [
    ("RNA_up", "RNA_log2FC", "RNA_pvalue_adjusted", 1),
    ("RNA_down", "RNA_log2FC", "RNA_pvalue_adjusted", -1),
    ("RIBO_up", "RIBO_log2FC", "RIBO_pvalue_adjusted", 1),
    ("RIBO_down", "RIBO_log2FC", "RIBO_pvalue_adjusted", -1),
    ("TE_up", "TE_log2FC", "TE_pvalue_adjusted", 1),
    ("TE_down", "TE_log2FC", "TE_pvalue_adjusted", -1),
]


def read_deltate_table(path, columns):
    """One deltaTE table, keyed by feature identifier.

    The files are tab separated with row names and a header one field shorter
    than the data rows, which is what makes pandas take the first column as the
    index. An empty file means deltaTE produced nothing, usually for lack of
    replicates, and is handled as an empty table rather than an error.
    """
    try:
        frame = pd.read_csv(path, sep="\t", comment="#")
        frame.index.name = "Identifier"
        return frame.reset_index(level=["Identifier"])
    except pd.errors.EmptyDataError:
        print("Warning deltaTE output is empty. Likely this is due to lacking replicates.")
        return pd.DataFrame(columns=["Identifier"] + columns)


def create_combined_dict(ribo_df, rna_df, te_df):
    """{identifier: [ribo row, rna row, te row]}, with None where a table lacks it."""
    combined = {}
    for position, frame in enumerate((ribo_df, rna_df, te_df)):
        for row in frame.itertuples(index=False):
            identifier = getattr(row, "Identifier")
            combined.setdefault(identifier, [None, None, None])[position] = row
    return combined


def deltate_output(args):
    genome_dict = eu.read_genome_dict(args.genome)
    annotation_dict = eu.annotation_to_dict(args.annotation_file)

    ribo_df = read_deltate_table(args.input_ribo, RIBO_RNA_COLUMNS)
    rna_df = read_deltate_table(args.input_rna, RIBO_RNA_COLUMNS)
    te_df = read_deltate_table(args.input_te, TE_COLUMNS)
    combined_dict = create_combined_dict(ribo_df, rna_df, te_df)

    header = (
        eu.DIFFEX_IDENTITY_HEADER
        + [header_name("RIBO", column) for column in RIBO_RNA_COLUMNS]
        + [header_name("RNA", column) for column in RIBO_RNA_COLUMNS]
        + [header_name("TE", column) for column in TE_COLUMNS]
        + ["Length", "Codon_count", "Start_codon", "Stop_codon", "Nucleotide_seq", "Aminoacid_seq"]
    )

    records = []
    for unique_id, rows in combined_dict.items():
        # A feature missing from any of the three tables cannot be reported.
        if None in rows:
            continue

        chromosome, start, stop, strand, gene_name, locus_tag, old_locus_tag = eu.resolve_location(
            unique_id, annotation_dict
        )

        start = int(start)
        stop = int(stop)
        length = stop - start + 1
        codon_count = int(length / 3)

        start_codon, stop_codon, nucleotide_seq, aa_seq = "", "", "", ""
        if chromosome in genome_dict:
            start_codon, stop_codon, nucleotide_seq, aa_seq, _ = eu.get_genome_information(
                genome_dict[chromosome], start - 1, stop - 1, strand
            )

        # Drop the leading Identifier field from each table's row.
        statistics = [value for row in rows for value in list(row)[1:]]

        records.append(
            [unique_id, chromosome, start, stop, strand, locus_tag, old_locus_tag, gene_name]
            + statistics
            + [length, codon_count, start_codon, stop_codon, nucleotide_seq, aa_seq]
        )

    all_df = pd.DataFrame.from_records(records, columns=header)
    all_df = all_df.sort_values(by=["TE_pvalue_adjusted", "Genome", "Start", "Stop", "Strand"])

    dataframe_dict = {"all": all_df}
    for sheet, log2fc_column, padj_column, direction in SPLITS:
        significant = all_df[padj_column] <= args.padj_cutoff
        if direction > 0:
            selected = all_df[log2fc_column] >= args.log2fc_cutoff
        else:
            selected = all_df[log2fc_column] <= args.log2fc_cutoff * -1
        dataframe_dict[sheet] = all_df[selected & significant]

    eu.excel_writer(args.output, dataframe_dict, [])


def main():
    parser = argparse.ArgumentParser(description="create excel files from deltaTE output")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_file",
                        required=True, help="annotation file")
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True,
                        help="reference genome")
    parser.add_argument("-i", "--delta_ribo", action="store", dest="input_ribo", required=True,
                        help="input txt file for ribo")
    parser.add_argument("-r", "--delta_rna", action="store", dest="input_rna", required=True,
                        help="input txt file for rna")
    parser.add_argument("-t", "--delta_te", action="store", dest="input_te", required=True,
                        help="input txt file for te")
    parser.add_argument("--padj_cutoff", action="store", dest="padj_cutoff", default=0.05,
                        type=float, help="padj cutoff for the differential expression analysis")
    parser.add_argument("--log2fc_cutoff", action="store", dest="log2fc_cutoff", default=1.0,
                        type=float, help="log2fc cutoff for the differential expression analysis")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True,
                        help="output xlsx file")
    args = parser.parse_args()

    deltate_output(args)


if __name__ == '__main__':
    main()
