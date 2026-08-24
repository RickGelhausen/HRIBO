#!/usr/bin/env python
"""Turn the riborex differential expression results into a spreadsheet.

Only the statistics columns and the cutoff columns are specific to riborex; the
rest is shared with the xtail and deltaTE tables in excel_utils.
"""

import argparse

import pandas as pd

import excel_utils as eu

# riborex writes DESeq2 results, so the columns are DESeq2's.
STATISTICS = [
    ("baseMean", "baseMean"),
    ("log2FoldChange", "log2FoldChange"),
    ("lfcSE", "lfcSE"),
    ("stat", "stat"),
    ("pvalue", "pvalue"),
    ("padj", "padj"),
]

# R's write.csv emits row names as an unnamed first column, which pandas reaches
# as _0.
IDENTIFIER_FIELD = "_0"


def riborex_output(args):
    genome_dict = eu.read_genome_dict(args.genome)
    annotation_dict = eu.annotation_to_dict(args.annotation_file)
    diff_expr_df = pd.read_csv(args.input_csv, sep=",", comment="#")

    all_df = eu.build_diffex_table(
        diff_expr_df, annotation_dict, genome_dict, STATISTICS, IDENTIFIER_FIELD
    )
    dataframe_dict = eu.split_up_down(
        all_df, "log2FoldChange", "padj", args.log2fc_cutoff, args.padj_cutoff
    )

    eu.excel_writer(args.output, dataframe_dict, [])


def main():
    parser = argparse.ArgumentParser(description="create excel files from riborex output")
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

    riborex_output(args)


if __name__ == '__main__':
    main()
