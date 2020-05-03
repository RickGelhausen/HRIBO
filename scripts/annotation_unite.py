#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections
import csv

def unite_annotation(args):
    """
    write the final annotation as gff
    """
    annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)
    nTuple = collections.namedtuple('Pandas', ["seq_name","source","feature","start","stop","score","strand","phase","attribute"])
    rows_unite = []
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        source = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        score = getattr(row, "_5")
        strand = getattr(row, "_6")
        phase = getattr(row, "_7")
        attribute = getattr(row, "_8")

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attribute) if x != ""]
        id = ""
        if "ID" in attribute_list:
            id = attribute_list[attribute_list.index("ID")+1]
        else:
            print("missing ID! Check your annotation.")
        attribute = "".join(["%s=%s;" % (attribute_list[i], attribute_list[i+1]) for i in range(0, len(attribute_list), 2)])

        rows_unite.append(nTuple(reference_name, source, feature, start, stop, score, strand, phase, attribute))

    return rows_unite


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='ensure that every row is in gff format.')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "input annotation file")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "output annotation file")
    args = parser.parse_args()
    rows_unite = unite_annotation(args)
    unite_df = pd.DataFrame.from_records(rows_unite, columns=[0,1,2,3,4,5,6,7,8])
    unite_df.to_csv(args.output, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
