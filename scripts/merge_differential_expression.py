#!/usr/bin/env python
'''This script takes input gff3 files and handles
overlapping intervals, by removing duplicates and
finding the longest non-overlapping interval.
'''

import pandas as pd
import argparse
import os, sys
import csv
import collections

def pool_contrasts_xtail(args):
    """
    collect all contrasts and pool files together
    """
    contrast_dict = {}
    for contrast_file in args.contrast_xlsx:
        cur_contrast = os.path.basename(contrast_file).split("_")[0]

        contrast_df = pd.read_excel(contrast_file, engine="openpyxl", sheet_name=0)

        for row in contrast_df.itertuples(index=False, name='Pandas'):
            gene_id = getattr(row, "Identifier")
            mRNA_log2FC = getattr(row, "mRNA_log2FC")
            RPF_log2FC = getattr(row, "RPF_log2FC")
            log2FC_te_v1 = getattr(row, "log2FC_TE_v1")
            pvalue_v1 = getattr(row, "pvalue_v1")
            log2FC_te_v2 = getattr(row, "log2FC_TE_v2")
            pvalue_v2 = getattr(row, "pvalue_v2")
            log2FC_te_final = getattr(row, "log2FC_TE_final")
            pvalue_final = getattr(row, "pvalue_final")
            pvalue_adj = getattr(row, "pvalue_adjusted")

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
    for contrast_file in args.contrast_xlsx:
        cur_contrast = os.path.basename(contrast_file).split("_")[0]

        contrast_df = pd.read_excel(contrast_file, engine="openpyxl", sheet_name=0)

        for row in contrast_df.itertuples(index=False, name='Pandas'):
            gene_id = getattr(row, "Identifier")
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

def pool_contrasts_deltate(args):
    """
    collect all contrasts and pool files together
    """
    contrast_dict = {}
    for contrast_file in args.contrast_xlsx:
        cur_contrast = os.path.basename(contrast_file).split("_")[0]

        contrast_df = pd.read_excel(contrast_file, engine="openpyxl", sheet_name=0)

        for row in contrast_df.itertuples(index=False, name='Pandas'):
            gene_id = getattr(row, "Identifier")

            ribo_baseMean = getattr(row, "RIBO_baseMean")
            ribo_log2FoldChange = getattr(row, "RIBO_log2FoldChange")
            ribo_lfcSE = getattr(row, "RIBO_lfcSE")
            ribo_pvalue = getattr(row, "RIBO_pvalue")
            ribo_padj = getattr(row, "RIBO_padj")

            rna_baseMean = getattr(row, "RNA_baseMean")
            rna_log2FoldChange = getattr(row, "RNA_log2FoldChange")
            rna_lfcSE = getattr(row, "RNA_lfcSE")
            rna_pvalue = getattr(row, "RNA_pvalue")
            rna_padj = getattr(row, "RNA_padj")

            te_baseMean = getattr(row, "TE_baseMean")
            te_log2FoldChange = getattr(row, "TE_log2FoldChange")
            te_lfcSE = getattr(row, "TE_lfcSE")
            te_stat = getattr(row, "TE_stat")
            te_pvalue = getattr(row, "TE_pvalue")
            te_padj = getattr(row, "TE_padj")

            contrast_tuple = (ribo_baseMean, ribo_log2FoldChange, ribo_lfcSE, ribo_pvalue, ribo_padj, \
                              rna_baseMean, rna_log2FoldChange, rna_lfcSE, rna_pvalue, rna_padj, \
                              te_baseMean, te_log2FoldChange, te_lfcSE, te_stat, te_pvalue, te_padj, \
                              args.tool + "_" + cur_contrast)

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

def merge_deltate_results(contrast_dict):
    """
    use the contrast dictionary to filter for the best result and add contrasts
    """

    header = ["gene_id", "RIBO_baseMean", "RIBO_log2FoldChange", "RIBO_lfcSE", "RIBO_pvalue", "RIBO_padj", \
              "RNA_baseMean", "RNA_log2FoldChange", "RNA_lfcSE", "RNA_pvalue", "RNA_padj", \
              "TE_baseMean", "TE_log2FoldChange", "TE_lfcSE", "TE_stat", "TE_pvalue", "TE_padj"]

    nTuple = collections.namedtuple('Pandas', header)

    rows = []
    for key, value in contrast_dict.items():
        for it in value:
            rows.append(nTuple(key, *it))

    return pd.DataFrame.from_records(rows, columns=header)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='pool the different contrasts together for each tool.')
    parser.add_argument("contrast_xlsx", nargs="*", metavar="contrast_xlsx", help="Path to input csv files.")
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
    elif args.tool == "deltate":
        deltate_dict = pool_contrasts_deltate(args)
        result_df = merge_deltate_results(deltate_dict)
    else:
        sys.exit("Invalid toolname!!!")

    with open(args.output_csv, "w") as f:
        result_df.to_csv(f, sep=",", index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
