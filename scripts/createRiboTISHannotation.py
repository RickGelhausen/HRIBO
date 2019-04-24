#!/usr/bin/env python

import os, sys
import argparse
import re
import pandas as pd
import csv
import collections

# helper function to create a dictionary {geneID : namedtuple}
def create_dictionary(input_df):
    print("Reading table rows!")

    gene_dict = dict()
    for row in input_df.itertuples(index=False, name='Pandas'):
        attributes = getattr(row, "_8")

        id = ""
        if "gene_id" in attributes.lower():
            attributes =  re.split('[; ]', attributes)
            id = attributes[attributes.index("gene_id")+1].replace("\"", "")

        elif "transcript_id" in attributes.lower():
            attributes =  re.split('[; ]', attributes)
            id = attributes[attributes.index("transcript_id")+1].replace("\"", "")

        if id != "":
            # save the feature and according row
            if id not in gene_dict:
                gene_dict[id] = (getattr(row,"_0"), getattr(row,"_3"), getattr(row,"_4"), getattr(row,"_5"),\
                                 getattr(row,"_6"), getattr(row,"_7"), 'gene_id "%s"; transcript_id "%s";' % (id, id))
    return gene_dict

# update annotation according to what is missing
def change_annotation(args):
    try:
        annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)
    except pd.errors.ParserError as error:
        sys.exit("Error reading the annotation file, please ensure that it is formatted correctly.\n"\
                +"Exiting with:\n   %s" % error)

    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])

    wanted_features = ["gene", "exon", "CDS"]
    gene_dict = create_dictionary(annotation_df)
    rows = []
    for key, val in gene_dict.items():
        for feature in wanted_features:
            rows.append(nTuple(val[0], "generated", feature, val[1], val[2], val[3], val[4], val[5], val[6]))

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create a dummy annotation file to make ribotish work')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True
                        , help="the annotation file to be validated.")
    parser.add_argument("--annotation_output", action="store", dest="annotation_output", required=True
                        , help="output annotation file.")

    args = parser.parse_args()

    df = change_annotation(args)

    df.to_csv(args.annotation_output, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
