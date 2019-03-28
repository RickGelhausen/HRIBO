#!/usr/bin/env python

import os, sys
import argparse
import re
import pandas as pd
import csv
import collections

# helper function to create a dictionary {geneID : namedtuple}
def create_dictionary(cleansedDF):
    print("Reading table rows!")
    geneDict = dict()
    insensitive_geneid= re.compile(re.escape("gene_id"), re.IGNORECASE)
    for row in cleansedDF.itertuples(index=False, name='Pandas'):
        commentPart = getattr(row, "_8")
        if "gene_id" in commentPart.lower():
            commentPart = commentPart.split(";")
            for opt in commentPart:
                if "gene_id" in opt.lower():
                    geneID = re.split('[ "]', opt)[2]

                    # save the row into the dictionary and ensure gene_id is written in lowercase
                    if geneID in geneDict:
                        geneDict[geneID].append(createNTuple(row, s8=insensitive_geneid.sub("gene_id", getattr(row, "_8"))))
                    else:
                        geneDict[geneID] = [createNTuple(row, s8=insensitive_geneid.sub("gene_id", getattr(row, "_8")))]
    return geneDict

# create a dataframe from the dictionary
def createAnnotationFromDictionary(geneDict):
    print("Creating new annotation file!")
    rows = []
    keys = sorted(geneDict.keys())
    for key in keys:
        for entry in sorted(geneDict[key], key=lambda x: x.s2):
            rows.append(entry)

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

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

# calculate the start codon for a given gene row
def generateStartCodon(row, newAttr):
    if getattr(row, "_6") == "+":
        newRightBoundary = str(int(getattr(row, "_3")) + 2)
        return createNTuple(row, s1="generated", s2="start_codon", s4=newRightBoundary,s8=newAttr)
    else:
        newLeftBoundary = str(int(getattr(row, "_4") - 2))
        return createNTuple(row, s1="generated", s2="start_codon", s3=newLeftBoundary,s8=newAttr)


# calculate the stop codon for a given cds row (ensembl style)
def generateStopCodon(row, newAttr):
    if getattr(row, "_6") == "+":
        newLeftBoundary = str(int(getattr(row, "_4") - 2))
        return createNTuple(row, s1="generated", s2="stop_codon", s3=newLeftBoundary,s8=newAttr)
    else:
        newRightBoundary = str(int(getattr(row, "_3")) + 2)
        return createNTuple(row, s1="generated", s2="stop_codon",s4=newRightBoundary,s8=newAttr)

# create exon
def generateExon(row, newAttr):
    return createNTuple(row, s1="generated", s2="exon", s8=newAttr)


def change_annotation(args):
    try:
        annotationDF = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)
    except pd.errors.ParserError as error:
        sys.exit("Error reading the annotation file, please ensure that it is formatted correctly.\n"\
                +"Exiting with:\n   %s" % error)

    rows = []
    for row in annotationDF.itertuples(index=False, name='Pandas'):
        attributes = getattr(row, "_8")
        if not "gene_id" in attributes:
            attrSplit = attributes.split(";")
            attributes = "gene_id " + attrSplit[0].split(" ")[1]+ "; " + attrSplit[0]

        rows.append(createNTuple(row, s2="exon", s8=attributes))
        rows.append(createNTuple(row, s8=attributes))

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create a dummy annotation file to make ribotish work')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True
                        , help="the annotation file to be validated.")
    parser.add_argument("--annotationOutput", action="store", dest="annotationOutput", required=True
                        , help="output annotation file.")

    args = parser.parse_args()

    df = change_annotation(args)

    df.to_csv(args.annotationOutput, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
