#!/usr/bin/env python
'''This script takes input gff3 files and handles
overlapping intervals, by removing duplicates and
finding the longest non-overlapping interval.
'''
from operator import itemgetter
import pandas as pd
import re
import argparse
import numpy as np
import os
import csv
import collections

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
            if ";" in attributes and "=" in attributes:
                attribute_list = [x for x in re.split('[;=]', attributes) if x != ""]
            else:
                attribute_list = [x.replace("\"", "") for x in re.split('[; ]', attributes) if x != ""]

            if len(attribute_list) % 2 == 0:
                for i in range(len(attribute_list)):
                    if i % 2 == 0:
                        attribute_list[i] = attribute_list[i].lower()
            else:
                print(attributes)
                sys.exit("Attributes section of gtf/gff is wrongly formatted!")

            key = ""
            if "id" in attribute_list:
                key = attribute_list[attribute_list.index("id")+1]

            gene = ""
            if "gene" in attribute_list:
                gene = attribute_list[attribute_list.index("gene")+1]

            locus_tag = "na"
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            gene_dict[key] = (gene, locus_tag)
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

        if ";" in attributes and "=" in attributes:
            attribute_list = [x for x in re.split('[;=]', attributes) if x != ""]
        else:
            attribute_list = [x.replace("\"", "") for x in re.split('[; ]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attributes)
            sys.exit("Attributes section of gtf/gff is wrongly formatted!")

        key = "%s:%s-%s:%s" % (reference_name, start, stop, strand)

        parent = ""
        if "parent" in attribute_list:
            parent = attribute_list[attribute_list.index("parent")+1]

        if parent in parent_dict:
            name, locus_tag = parent_dict[parent]

            if name == "":
                name = "%s:%s-%s:%s" % (reference_name, start, stop, strand)
                if "name" in attribute_list:
                    name = attribute_list[attribute_list.index("name")+1]
                if "gene" in attribute_list:
                    name = attribute_list[attribute_list.index("gene")+1]

            if locus_tag == "":
                locus_tag = "na"
                if "locus_tag" in attribute_list:
                    locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

        else:
            name = "%s:%s-%s:%s" % (reference_name, start, stop, strand)
            if "name" in attribute_list:
                name = attribute_list[attribute_list.index("name")+1]
            if "gene" in attribute_list:
                name = attribute_list[attribute_list.index("gene")+1]

            locus_tag = "na"
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

        annotation_dict[key] = (name, locus_tag)

    return annotation_dict

def generate_output_gff(args, overlap_dict):
    """
    write an output file where only the best rank is taken for each overlapping prediction
    """
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])
    annotation_dict = generate_annotation_dict(args)

    rows = []
    for key, value in overlap_dict.items():
        reference_name, mid, strand = key.split(":")
        start, stop = mid.split("-")
        evidence = set()
        cur_rank = 999999
        cur_pred_value = -10000.0
        for pred, dist, attribute in value:
            if ";" in attribute and "=" in attribute:
                attribute_list = [x for x in re.split('[;=]', attribute) if x != ""]
            else:
                attribute_list = [x.replace("\"", "") for x in re.split('[; ]', attribute) if x != ""]

            if len(attribute_list) % 2 == 0:
                for i in range(len(attribute_list)):
                    if i % 2 == 0:
                        attribute_list[i] = attribute_list[i].lower()
            else:
                print(attribute)
                sys.exit("Attributes section of gtf/gff is wrongly formatted!")

            pred_value = pred
            if pred_value >= cur_pred_value:
                cur_pred_value = pred_value

            if "condition" in attribute_list and "method" in attribute_list and "replicate" in attribute_list:
                condition = attribute_list[attribute_list.index("condition")+1]
                #method = attribute_list[attribute_list.index("method")+1]
                replicate = attribute_list[attribute_list.index("replicate")+1]
                evidence.add(condition + "-" + replicate)
            elif "condition" in attribute_list and "method" in attribute_list:
                condition = attribute_list[attribute_list.index("condition")+1]
                method = attribute_list[attribute_list.index("method")+1]
                evidence.add(method + "-" + condition)

            if key in annotation_dict:
                name, locus_tag = annotation_dict[key]
            else:
                name, locus_tag = key, "na"

            new_attributes = "ID=%s;Name=%s;Locus_tag=%s;Pred_value=%s;Evidence=%s;" % (key, name,locus_tag,cur_pred_value, " ".join(evidence))
 
        rows.append(nTuple(reference_name, "merged", "CDS", start, stop, cur_pred_value, strand, dist, new_attributes))
    return  pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])

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
        newDF = generate_output_gff(args, orf_dict)
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


if __name__ == '__main__':
    main()
