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

# little helper function to create named tuple without having to always state every argument
def createNTuple(row, s0=None, s1=None, s2=None, s3=None, s4=None, s5=None, s6=None, s7=None, s8=None):
    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])
    try:
        c0, c1, c2 = getattr(row, "_0"), getattr(row, "_1"), getattr(row, "_2")
        c3, c4, c5 = getattr(row, "_3"), getattr(row, "_4"), getattr(row, "_5")
        c6, c7, c8 = getattr(row, "_6"), getattr(row, "_7"), getattr(row, "_8")
    except AttributeError:
        c0, c1, c2 = getattr(row, "s0"), getattr(row, "s1"), getattr(row, "s2")
        c3, c4, c5 = getattr(row, "s3"), getattr(row, "s4"), getattr(row, "s5")
        c6, c7, c8 = getattr(row, "s6"), getattr(row, "s7"), getattr(row, "s8")

    if s0 != None:
        c0 = s0
    if s1 != None:
        c1 = s1
    if s2 != None:
        c2 = s2
    if s3 != None:
        c3 = s3
    if s4 != None:
        c4 = s4
    if s5 != None:
        c5 = s5
    if s6 != None:
        c6 = s6
    if s7 != None:
        c7 = s7
    if s8 != None:
        c8 = s8
    return nTuple(c0, c1, c2, c3, c4, c5, c6, c7, c8)


# helper function to create a dictionary {geneID : namedtuple}
def create_dictionary(inputDF):
    geneDict = dict()
    for row in inputDF.itertuples(index=False, name='Pandas'):
        attributes = re.split('[;=]', getattr(row, "_8"))
        if "ID" in attributes:
            geneID = attributes[attributes.index("ID")+1]

            # save the row into the dictionary and ensure gene_id is written in lowercase
            if geneID in geneDict:
                geneDict[geneID].append(createNTuple(row))
            else:
                geneDict[geneID] = [createNTuple(row)]
    return geneDict

def handle_overlap(args):
    inputDF = pd.read_csv(args.inputGFF, sep='\t', header=None)

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


        attribute = "ID="+key+";Name="+key+";ORF_type="+",".join(orftype)+";Evidence="+",".join(evidence)
        rows.append(createNTuple(sampleRow,s1="merged",s8=attribute))

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
        with open(args.outputGFF, "w") as f:
            f.write("##gff-version 3\n")
            newDF = handle_overlap(args)
            newDF.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)




if __name__ == '__main__':
    main()
