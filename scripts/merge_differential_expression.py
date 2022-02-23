#!/usr/bin/env python
'''This script takes input gff3 files and handles
overlapping intervals, by removing duplicates and
finding the longest non-overlapping interval.
'''

import pandas as pd
import argparse
import os
import csv
import collections

def pool_contrasts_xtail(args):
    """
    collect all contrasts and pool files together
    """
    contrast_dict = {}
    for contrast_file in args.contrast_csv:
        cur_contrast = os.path.basename(contrast_file).split("_")[0]

        contrast_df = pd.read_csv(contrast_file, sep=",", comment="#")

        for row in contrast_df.itertuples(index=False, name='Pandas'):
            gene_id = getattr(row, "_0")
            mRNA_log2FC = getattr(row, "mRNA_log2FC")
            RPF_log2FC = getattr(row, "RPF_log2FC")
            log2FC_te_v1 = getattr(row, "log2FC_TE_v1")
            pvalue_v1 = getattr(row, "pvalue_v1")
            log2FC_te_v2 = getattr(row, "log2FC_TE_v2")
            pvalue_v2 = getattr(row, "pvalue_v2")
            log2FC_te_final = getattr(row, "log2FC_TE_final")
            pvalue_final = getattr(row, "pvalue_final")
            pvalue_adj = getattr(row, "_9")

            contrast_tuple = (mRNA_log2FC, RPF_log2FC, log2FC_te_v1, pvalue_v1, log2FC_te_v2, pvalue_v2, \
                              log2FC_te_final, pvalue_final, pvalue_adj, args.tool + "_" + cur_contrast)

            if gene_id in contrast_dict:
                contrast_dict[gene_id].append(contrast_tuple)
            else:
                contrast_dict[gene_id] = [contrast_tuple]

    return contrast_dict

def pool_contrasts_riborex(args):
    """
    collect all contrasts and pool files together
    """
    contrast_dict = {}
    for contrast_file in args.contrast_csv:
        cur_contrast = os.path.basename(contrast_file).split("_")[0]

        contrast_df = pd.read_csv(contrast_file, sep=",", comment="#")

        for row in contrast_df.itertuples(index=False, name='Pandas'):
            gene_id = getattr(row, "_0")
            baseMean = getattr(row, "baseMean")
            log2FoldChange = getattr(row, "log2FoldChange")
            lfcSE = getattr(row, "lfcSE")
            stat = getattr(row, "stat")
            pvalue = getattr(row, "pvalue")
            padj = getattr(row, "padj")
            contrast_tuple = (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, args.tool + "_" + cur_contrast)

            if gene_id in contrast_dict:
                contrast_dict[gene_id].append(contrast_tuple)
            else:
                contrast_dict[gene_id] = [contrast_tuple]

    return contrast_dict

def merge_xtail_results(contrast_dict):
    """
    use the contrast dictionary to filter for the best result and add contrasts
    """

    header = ["gene_id", "mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1", \
              "log2FC_TE_v2", "pvalue_v2", "log2FC_TE_final", "pvalue_final", \
              "pvalue_adjust", "contrast"]

    nTuple = collections.namedtuple('Pandas', header)

    rows = []
    for key, value in contrast_dict.items():
        for it in value:
            rows.append(nTuple(key, *it))

    return pd.DataFrame.from_records(rows, columns=header)

def merge_riborex_results(contrast_dict):
    """
    use the contrast dictionary to filter for the best result and add contrasts
    """

    header = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "contrast"]

    nTuple = collections.namedtuple('Pandas', header)

    rows = []
    for key, value in contrast_dict.items():
        for it in value:
            rows.append(nTuple(key, *it))

    return pd.DataFrame.from_records(rows, columns=header)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='pool the different contrasts together for each tool.')
    parser.add_argument("contrast_csv", nargs="*", metavar="contrast_csv", help="Path to input csv files.")
    parser.add_argument("-o", "--output_csv", action="store", dest="output_csv", required=True
                                           , help= "the output csv file.")
    parser.add_argument("-t", "--tool", action="store", dest="tool", required=True, help= "tool: riborex/xtail")
    args = parser.parse_args()

    if args.tool == "riborex":
        riborex_dict = pool_contrasts_riborex(args)
        result_df = merge_riborex_results(riborex_dict)
    elif args.tool == "xtail":
        xtail_dict = pool_contrasts_xtail(args)
        result_df = merge_xtail_results(xtail_dict)
    else:
        sys.exit("invalid toolname!!!")

    with open(args.output_csv, "w") as f:
        result_df.to_csv(f, sep=",", index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
