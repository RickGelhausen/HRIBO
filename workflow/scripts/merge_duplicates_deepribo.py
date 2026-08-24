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
import sys

import gff_utils

def generate_dictionary(args):
    """
    read the input file and create a dictionary containing information on overlapping genes
    (rank, attributes)
    """
    in_df = pd.read_csv(args.inputGFF, comment="#", sep="\t", header=None)

    overlap_dict = {}
    for row in in_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        source = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        prediction_rank = getattr(row, "_5")
        strand = getattr(row, "_6")
        dist = getattr(row, "_7")
        attributes = getattr(row, "_8")

        key = "%s:%s-%s:%s" % (reference_name, start, stop, strand)
        if key in overlap_dict:
            overlap_dict[key].append((prediction_rank, dist, attributes))
        else:
            overlap_dict[key] = [(prediction_rank, dist, attributes)]

    return overlap_dict

def create_gene_dict(annotation_df):
    gene_dict = {}
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        feature = getattr(row, "_2")
        attributes = getattr(row, "_8")
        if feature == "gene":
            parsed = gff_utils.parse_attributes(attributes)

            gene_dict[gff_utils.first_attribute(parsed, "id")] = (
                gff_utils.first_attribute(parsed, "gene"),
                gff_utils.first_attribute(parsed, "locus_tag"),
                gff_utils.first_attribute(parsed, "old_locus_tag"),
            )
    return gene_dict


def generate_annotation_dict(args):
    annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    parent_dict = create_gene_dict(annotation_df)

    annotation_dict = {}
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")

        if feature not in ["CDS", "cds"]:
            continue

        parsed = gff_utils.parse_attributes(attributes)

        key = "%s:%s-%s:%s" % (reference_name, start, stop, strand)

        parent = gff_utils.first_attribute(parsed, "parent")

        # The feature's own attributes, used where the gene feature has none.
        own_name = gff_utils.first_attribute(parsed, "name", "gene", default=key)
        own_locus_tag = gff_utils.first_attribute(parsed, "locus_tag")
        # Previously read the "locus_tag" attribute here, so old_locus_tag was
        # filled with the current locus tag instead of the old one.
        own_old_locus_tag = gff_utils.first_attribute(parsed, "old_locus_tag")

        if parent in parent_dict:
            name, locus_tag, old_locus_tag = parent_dict[parent]
            if name == "":
                name = own_name
            if locus_tag == "":
                locus_tag = own_locus_tag or "na"
            if old_locus_tag == "":
                old_locus_tag = own_old_locus_tag
        else:
            # old_locus_tag was never initialised on this branch, so it either
            # raised or silently carried over the previous row's value.
            name = own_name
            locus_tag = own_locus_tag
            old_locus_tag = own_old_locus_tag

        annotation_dict[key] = (name, locus_tag, old_locus_tag)

    return annotation_dict

def generate_output_gff(args, overlap_dict):
    """
    write an output file where only the best rank is taken for each overlapping prediction
    """
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])
    annotation_dict = generate_annotation_dict(args)

    rows = []
    rows_plus = []
    for key, value in overlap_dict.items():
        reference_name, mid, strand = key.split(":")
        start, stop = mid.split("-")
        evidence = set()
        cur_rank = 999999
        cur_pred_value = -10000.0
        for pred, dist, attribute in value:
            # Handles both the GFF3 and the GTF2 attribute forms.
            parsed = gff_utils.parse_attributes(attribute)

            pred_value = pred
            if pred_value >= cur_pred_value:
                cur_pred_value = pred_value

            # A replicate identifies the evidence precisely; without one, the
            # method and condition are the best that can be said.
            if {"condition", "method", "replicate"} <= parsed.keys():
                evidence.add(parsed["condition"] + "-" + parsed["replicate"])
            elif {"condition", "method"} <= parsed.keys():
                evidence.add(parsed["method"] + "-" + parsed["condition"])

            if key in annotation_dict:
                name, locus_tag, old_locus_tag = annotation_dict[key]
            else:
                name, locus_tag, old_locus_tag = key, "", ""


            new_attributes = "ID=%s;Name=%s;" % (key, name)
            if locus_tag != "":
                new_attributes += "locus_tag=%s;" % (locus_tag)
            if old_locus_tag != "":
                new_attributes += "old_locus_tag=%s;" % (old_locus_tag)

            # sorted(): see merge_duplicates_reparation.py; an unsorted set join makes
            # the output vary between runs.
            new_attributes += "Pred_value=%s;Evidence=%s;" % (cur_pred_value, " ".join(sorted(evidence)))

        if cur_pred_value >= 0:
            rows_plus.append(nTuple(reference_name, "deepribo", "CDS", start, stop, cur_pred_value, strand, dist, new_attributes))
        rows.append(nTuple(reference_name, "deepribo", "CDS", start, stop, cur_pred_value, strand, dist, new_attributes))

    return pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"]), \
           pd.DataFrame.from_records(rows_plus, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='condense duplicates into one entry')
    parser.add_argument("-i", "--inputGFF", action="store", dest="inputGFF", required=True
                                          , help= "the input file (gff3 format).")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True
                                           , help= "annotation file")
    args = parser.parse_args()

    if os.stat(args.inputGFF).st_size == 0:
        open(args.outputGFF, 'a').close()
    else:
        orf_dict = generate_dictionary(args)
        newDF, plusDF = generate_output_gff(args, orf_dict)
        newDF = newDF.sort_values(by=["score"], ascending=False)
        newDF = newDF.reset_index()
        dist_list = list(newDF["phase"])

        counter = 1
        for i in range(len(dist_list)):
            if dist_list[i] == -1:  ## SOME ENTRIES ARE NEITHER 0 nor -1, they are not considered novel due to a lack of definition
                dist_list[i] = counter
                counter += 1
            else:
                dist_list[i] = 999999

        newDF["phase"] = dist_list

        newDF["score"] = newDF.index + 1
        newDF = newDF.drop(columns=["index"])

        with open(args.outputGFF, "w") as f:
            f.write("##gff-version 3\n")
        with open(args.outputGFF, "a") as f:
            newDF.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

        with open(args.outputGFF.replace(".gff", "_plus.gff"), "w") as f:
            f.write("##gff-version 3\n")
        with open(args.outputGFF.replace(".gff", "_plus.gff"), "a") as f:
            plusDF.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
