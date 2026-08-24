#!/usr/bin/env python
"""Turn the xtail differential expression results into a spreadsheet.

Only the statistics columns and the cutoff columns are specific to xtail; the
rest is shared with the riborex and deltaTE tables in excel_utils.
"""

import argparse

import pandas as pd

import excel_utils as eu

# xtail reports two variants of the TE fold change plus a combined final value.
# The last column is R-named "pvalue.adjust", which is not a valid Python
# identifier, so pandas exposes it positionally as _9.
STATISTICS = [
    ("mRNA_log2FC", "mRNA_log2FC"),
    ("RPF_log2FC", "RPF_log2FC"),
    ("log2FC_TE_v1", "log2FC_TE_v1"),
    ("pvalue_v1", "pvalue_v1"),
    ("log2FC_TE_v2", "log2FC_TE_v2"),
    ("pvalue_v2", "pvalue_v2"),
    ("log2FC_TE_final", "log2FC_TE_final"),
    ("pvalue_final", "pvalue_final"),
    ("pvalue_adjusted", "_9"),
]

# R's write.csv emits row names as an unnamed first column.
IDENTIFIER_FIELD = "_0"


def xtail_output(args):
    genome_dict = eu.read_genome_dict(args.genome)
    annotation_dict = eu.annotation_to_dict(args.annotation_file)
    diff_expr_df = pd.read_csv(args.input_csv, sep=",", comment="#")

    all_df = eu.build_diffex_table(
        diff_expr_df, annotation_dict, genome_dict, STATISTICS, IDENTIFIER_FIELD
    )
    dataframe_dict = eu.split_up_down(
        all_df, "log2FC_TE_final", "pvalue_adjusted", args.log2fc_cutoff, args.padj_cutoff
    )

    eu.excel_writer(args.output, dataframe_dict, [])


def main():
    parser = argparse.ArgumentParser(description="create excel files from xtail output")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_file",
                        required=True, help="annotation file")
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True,
                        help="reference genome")
    parser.add_argument("-i", "--input", action="store", dest="input_csv", required=True,
                        help="input csv file")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True,
                        help="output xlsx file")
    parser.add_argument("--padj_cutoff", action="store", dest="padj_cutoff", default=0.05,
                        type=float, help="padj cutoff for the differential expression analysis")
    parser.add_argument("--log2fc_cutoff", action="store", dest="log2fc_cutoff", default=1.0,
                        type=float, help="log2fc cutoff for the differential expression analysis")
    args = parser.parse_args()

    xtail_output(args)


if __name__ == '__main__':
    main()
