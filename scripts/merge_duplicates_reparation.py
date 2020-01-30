#!/usr/bin/env python
'''This script takes input gff3 files and handles
overlapping intervals, by merging duplicates.
'''
from operator import itemgetter
import pandas as pd
import re
import argparse
import numpy as np
import os
import csv
import collections


def create_dictionary(inputDF):
    """
    function to create a dictionary {geneID : namedtuple}
    """
    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])

    geneDict = dict()
    for row in inputDF.itertuples(index=False, name='Pandas'):
        attributes = re.split('[;=]', getattr(row, "_8"))
        if "ID" in attributes:
            geneID = attributes[attributes.index("ID")+1]
                # save the row into the dictionary and ensure gene_id is written in lowercase
            if geneID in geneDict:
                geneDict[geneID].append(nTuple(*row))
            else:
                geneDict[geneID] = [nTuple(*row)]

    return geneDict

def handle_overlap(args):
    """
    read the input gff and merge all duplicate intervals
    """
    inputDF = pd.read_csv(args.inputGFF, sep='\t', header=None)
    nTuple = collections.namedtuple('Pandas', ["seq_name","source","feature","start","stop","score","strand","phase","attribute"])

    # create a dictionary for common ids
    geneDict = create_dictionary(inputDF)

    # run over all entries in the dictionary and combine overlapping ones
    rows = []
    for key in geneDict.keys():
        sampleRow = geneDict[key][0]
        evidence = set()
        orftype = set()
        for row in geneDict[key]:
            attributes = re.split('[;=]', getattr(row, "s8"))
            if "Condition" in attributes and "Method" in attributes and "Replicate" in attributes:
                condition = attributes[attributes.index("Condition")+1]
                method = attributes[attributes.index("Method")+1]
                replicate = attributes[attributes.index("Replicate")+1]
                evidence.add(method + "-" + condition + "-" + replicate)
            elif "Condition" in attributes and "Method" in attributes:
                condition = attributes[attributes.index("Condition")+1]
                method = attributes[attributes.index("Method")+1]
                evidence.add(method + "-" + condition)

            if "ORF_type" in attributes:
                orftype.add(attributes[attributes.index("ORF_type")+1])


        attribute = "ID="+key+";Name="+key+";ORF_type="+",".join(orftype)+";Evidence="+" ".join(evidence)

        rows.append(nTuple(getattr(sampleRow, "s0"),"merged", getattr(sampleRow, "s2"), getattr(sampleRow, "s3"), \
                           getattr(sampleRow, "s4"), getattr(sampleRow, "s5"),getattr(sampleRow, "s6"), \
                           getattr(sampleRow, "s7"), attribute))

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='condense duplicates into one entry')
    parser.add_argument("-i", "--inputGFF", action="store", dest="inputGFF", required=True
                                          , help= "the input file (gff3 format).")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")
    args = parser.parse_args()
    if os.stat(args.inputGFF).st_size == 0:
        open(args.outputGFF, 'a').close()
    else:
        newDF = handle_overlap(args)
        newDF.to_csv(args.outputGFF, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)




if __name__ == '__main__':
    main()
