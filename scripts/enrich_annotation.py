#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import collections
import csv

def create_parent_dictionary(annotation_df):
    """
    collect all attributes of parents in a dictionary
    key: ID, value: attribute_list
    """
    parent_dict = {}
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        attributes = getattr(row, "_8")

        if "id=" not in attributes.lower():
            continue
        else:
            attribute_list = [x for x in re.split('[;=]', attributes)]
            try:
                idx = attribute_list[next(i for i,v in enumerate(attribute_list) if v.lower() == "id") + 1]
            except:
                print("Missing ID in row!")
                print(row)
                sys.exit()

            parent_dict[idx] = attribute_list

    return parent_dict


def enrich_children(annotation_df):
    """
    if no locus_tag or name is defined, get the one from the parent or the grand-parents ...
    """
    parent_dict = create_parent_dictionary(annotation_df)
    nTuple = collections.namedtuple('Pandas', ["seq_name","source","feature","start","stop","score","strand","phase","attributes"])

    new_rows = []
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        source = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        score = getattr(row, "_5")
        strand = getattr(row, "_6")
        phase = getattr(row, "_7")
        attributes = getattr(row, "_8")
        if "parent=" in attributes.lower():
            attribute_list = [x for x in re.split('[;=]', attributes)]

            has_parent = True
            while has_parent:
                try:
                    parent = attribute_list[next(i for i,v in enumerate(attribute_list) if v.lower() == "parent")+1]
                except StopIteration:
                    has_parent = False

                attribute_list = parent_dict[parent]


            if parent not in parent_dict:
                print("Warning! Missing parent!")
                print(row)
                new_rows.append(row)
                continue

            locus_tag = ""
            if "locus_tag" not in attributes.lower():
                try:
                    locus_tag = "locus_tag=%s;" % parent_dict[parent][next(i for i,v in enumerate(parent_dict[parent]) if v.lower() == "locus_tag")+1]
                except StopIteration:
                    locus_tag = ""

            old_locus_tag = ""
            if "old_locus_tag" not in attributes.lower():
                try:
                    old_locus_tag = "old_locus_tag=%s;" % parent_dict[parent][next(i for i,v in enumerate(parent_dict[parent]) if v.lower() == "old_locus_tag")+1]
                except StopIteration:
                    old_locus_tag = ""

            name = ""
            if "name" not in attributes.lower():
                try:
                    name = "Name=%s;" % parent_dict[parent][next(i for i,v in enumerate(parent_dict[parent]) if v.lower() == "name")+1]
                except:
                    name = ""

            new_rows.append(nTuple(reference_name, source, feature, start, stop, score, strand, phase, attributes + ";" + locus_tag + old_locus_tag + name))

        else:
            new_rows.append(row)

    return pd.DataFrame.from_records(new_rows, columns=[0,1,2,3,4,5,6,7,8])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='enrich the child features with the information provided in the parent.')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "input annotation (gff3 format)")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "output annotation, enriched with parent info")
    args = parser.parse_args()

    annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    with open(args.output, "w") as f:
        f.write("##gff-version 3\n")

    with open(args.output, "a") as f:
        annotation_df = enrich_children(annotation_df)
        annotation_df.to_csv(f, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
