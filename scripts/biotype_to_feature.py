#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import csv
import collections


def extract_biotype(args):
    """
    go through gtf2 file and add features for rRNA, tRNA, sRNA and ncRNA if detected in biotype
    """
    # read annotation file
    annotation_df = pd.read_csv(args.annotation, comment="#", sep="\t", header=None)
    annotation_dict = {}

    nTuple = collections.namedtuple('Pandas', ["seq_name","info","feature","start","stop","score","strand","phase","attributes"])

    rows = []
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        seq_name = str(getattr(row, "_0"))
        info = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = str(getattr(row, "_3"))
        stop = str(getattr(row, "_4"))
        score = str(getattr(row, "_5"))
        strand = str(getattr(row, "_6"))
        phase = str(getattr(row, "_7"))
        description = getattr(row, "_8")
        attributes = [x.replace("\"","") for x in re.split('[; ]', description) if x != ""]

        rows.append(row)
        if feature == "CDS" and "gene_biotype" in attributes:
            gene_biotype = attributes[attributes.index("gene_biotype") + 1]
            if gene_biotype in ["ncRNA", "sRNA", "tRNA", "rRNA", "mRNA"]:
                rows.append(nTuple(seq_name, "generated", gene_biotype, start, stop, score, strand, phase, description))

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='clean annotation gff3.')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "The standard annotation file.")
    parser.add_argument("-o", "--output", action="store", dest="outputGFF", required=True, help= "The cleaned output file.")
    args = parser.parse_args()

    newDF = extract_biotype(args)
    newDF.to_csv(args.outputGFF, header=None, sep="\t", index=False, mode="a", quoting=csv.QUOTE_NONE)



if __name__ == '__main__':
    main()
