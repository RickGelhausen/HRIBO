#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import csv
import collections

def parse_gff3(annotation_df):
    """
    parse the gff3 file into a list of information
    """

    entries = []
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        origin = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        score = getattr(row, "_5")
        strand = getattr(row, "_6")
        phase = getattr(row, "_7")
        attributes = [x for x in re.split('[;=]', getattr(row, "_8"))]

        if feature not in ["start_codon", "stop_codon"]:
            entries.append([reference_name, origin, feature, start, stop, score, strand, phase, attributes])

    return entries

def parse_gff2(annotation_df):
    """
    parse the gff2 file into a list of information
    """

    entries = []
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        origin = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        score = getattr(row, "_5")
        strand = getattr(row, "_6")
        phase = getattr(row, "_7")
        attributes = [x.replace("\"", "") for x in re.split('[; ]', getattr(row, "_8")) if x != ""]

        if feature not in ["start_codon", "stop_codon"]:
            entries.append([reference_name, origin, feature, start, stop, score, strand, phase, attributes])

    return entries


def create_annotation_dict(args):
    """
    create a dictionary out of the annotation file
    key : (reference, start, stop, strand)
    value: [unique_features,...]
    """

    annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    rows_with_id = sum(annotation_df[8].str.contains("ID=").tolist())
    number_of_rows = len(annotation_df.index)

    # gff3
    if rows_with_id >= number_of_rows * 0.75:
        annotation_rows = parse_gff3(annotation_df)
    else: #gff2
        annotation_rows = parse_gff2(annotation_df)

    annotation_dict = {}
    for row in annotation_rows:
        reference_name = row[0]
        origin = row[1]
        feature = row[2]
        start = row[3]
        stop = row[4]
        score = row[5]
        strand = row[6]
        phase = row[7]
        attributes = row[8]

        key = (reference_name, start, stop, strand)
        if key in annotation_dict:
            annotation_dict[key].append(feature)
        else:
            annotation_dict[key] = [feature]

    return annotation_dict

def generate_shortened_annotation(annotation_dict):
    nTuple = collections.namedtuple('Pandas', ["seq_name","info","feature","start","stop","score","strand","phase","attributes"])

    rows = []
    rrnaCounter=1
    trnaCounter=1
    ncrnaCounter=1
    geneCounter=1
    srnaCounter=1
    for key, features in annotation_dict.items():
        if "rRNA" in features:
            rows.append(nTuple(key[0],"extracted", "rRNA", key[1], key[2],".",key[3],".", "gene_id \"rrna%s\";" % rrnaCounter))
            rrnaCounter += 1
        elif "tRNA" in features:
            rows.append(nTuple(key[0],"extracted", "tRNA", key[1], key[2],".",key[3],".", "gene_id \"trna%s\";" % trnaCounter))
            trnaCounter += 1
        elif "sRNA" in features:
            rows.append(nTuple(key[0],"extracted", "sRNA", key[1], key[2],".",key[3],".", "gene_id \"srna%s\";" % srnaCounter))
            srnaCounter += 1
        elif "ncRNA" in features:
            rows.append(nTuple(key[0],"extracted", "ncRNA", key[1], key[2],".",key[3],".", "gene_id \"ncrna%s\";" % ncrnaCounter))
            ncrnaCounter += 1
        elif "gene" in features:
            rows.append(nTuple(key[0],"extracted", "gene", key[1], key[2],".",key[3],".", "gene_id \"gene%s\";" % geneCounter))
            geneCounter += 1
        else:
            print("Unexpected feature: ensure that gene features are available if no rRNA, tRNA, sRNA, ncRNA is referenced.")

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='extract only the rows interesting for featurecounts and generate a gtf2 file')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "the input annotation file.")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "the processed annotation file.")
    args = parser.parse_args()

    annotation_dict = create_annotation_dict(args)
    short_annotation_df = generate_shortened_annotation(annotation_dict)

    with open(args.output, "w") as f:
        f.write("##gff-version 2\n")
    with open(args.output, "a") as f:
        short_annotation_df.sort_values(by=[0,3,4], inplace=True)
        short_annotation_df.to_csv(f, header=None, sep="\t", index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
