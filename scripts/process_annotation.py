#!/usr/bin/env python
import csv
import pandas as pd
import argparse
import re
import collections
import fileinput

# Reformat the attribute column exactly how metagene wants it
def combineSplit(attrSplit):
    attribute = ""
    attrCount = 0
    openQuotes = False
    for attr in attrSplit:
        if openQuotes:
            if "\"" in attr:
                attribute += attr + "; "
                attrCount = 0
                openQuotes = False
            else:
                attribute += attr + " "

        else:
            if attrCount < 1:
                attrCount += 1
                attribute += attr + " "
            else:
                if attr.count("\"") == 2:
                    attribute += attr + "; "
                    attrCount = 0
                else:
                    attribute += attr + " "
                    openQuotes = True

    return attribute

def process_annotation(args):
    nTuple = collections.namedtuple('Pandas', ["seq_name","source","feature","start","stop","score","strand","phase","attribute"])
    allowedFeatures= ["exon", "cds", "start_codon", "stop_codon", "gene", "transcript"]
    # read input annotation
    annDF = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    # check if gff2/gtf or gff3 format

    # number of entries containing ID=
    trueValues = sum(annDF[8].str.contains("ID=").tolist())
    # number of rows
    nrows = len(annDF.index)
    # if 75% of all rows contain ID= it is likely a gff3
    if trueValues >= nrows * 0.75:
        # gff3
        with open(args.output, "w") as f:
            f.write( "##gff-version 3\n")
        with open(args.output, "a") as f:
            annDF.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    else:
        # gff2
        # only accept rows that contain gene_id or transcript_id, then generate the missing
        annDF = annDF[annDF[8].str.contains("gene_id") | annDF[8].str.contains("transcript_id")]

        # ensure that the separation of attributes is done correctly
        rows = []
        for row in annDF.itertuples(index=False, name='Pandas'):
            # ignore entries that do not contain wanted features
            if getattr(row, "_2").lower() not in allowedFeatures:
                continue

            attributes = getattr(row, "_8")
            attrSplit = list(filter(None, re.split('[ ;]', attributes)))

            # if either gene_id or transcript_id is missing, add an empty one
            if "gene_id" in attributes:
                if not "transcript_id" in attributes:
                    attrSplit = attrSplit[:2] + ['transcript_id','""'] + attrSplit[2:]
            elif "transcript_id" in attributes:
                if not "gene_id" in attributes:
                    attrSplit = ['gene_id','""'] + attrSplit

            attributes = combineSplit(attrSplit)
            rows.append(nTuple(getattr(row, "_0"), getattr(row, "_1"), getattr(row, "_2"), getattr(row, "_3"), \
                               getattr(row, "_4"), getattr(row, "_5"), getattr(row, "_6"), getattr(row, "_7"), \
                               attributes))

        annDF = pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

        # write processed annotation to output
        with open(args.output, "w") as f:
            f.write( "##gff-version 2\n")
        with open(args.output, "a") as f:
            annDF.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='process the annotation for use with plastid')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "the input annotation file.")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "the processed annotation file.")
    args = parser.parse_args()

    process_annotation(args)


if __name__ == '__main__':
    main()
