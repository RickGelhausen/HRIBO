#!/usr/bin/env python
'''This script takes an annotation file and a genome file as input
and validates them. If they are not valid, try creating valid ones.
(only keeping needed identifiers and creating start/stop codons like ensembl)
'''

import os, sys
import argparse
import re
import pandas as pd
import csv
import collections

# split the string into strings of given length
def splitString(s, l):
    return "".join([s[i:i+l] + "\n" for i in range(0, len(s), l)])

# calculate the start codon for a given gene row
def generateStartCodon(row):
    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])
    if getattr(row, "_6") == "+":
        newRightBoundary = str(int(getattr(row, "_3")) + 2)
        return nTuple(getattr(row, "_0"), getattr(row, "_1"), "start_codon"\
                            , getattr(row, "_3"), newRightBoundary\
                            , getattr(row, "_5"), getattr(row, "_6")\
                            , getattr(row, "_7"), getattr(row, "_8"))
    else:
        newLeftBoundary = str(int(getattr(row, "_4") - 2))
        return nTuple(getattr(row, "_0"), getattr(row, "_1"), "start_codon"\
                            , newLeftBoundary, getattr(row, "_4")\
                            , getattr(row, "_5"), getattr(row, "_6")\
                            , getattr(row, "_7"), getattr(row, "_8"))

# calculate the stop codon for a given cds row (ensembl style)
def generateStopCodon(row):
    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])
    if getattr(row, "_6") == "+":
        newLeftBoundary = str(int(getattr(row, "_4") - 2))
        return nTuple(getattr(row, "_0"), getattr(row, "_1"), "stop_codon"\
                            , newLeftBoundary, getattr(row, "_4")\
                            , getattr(row, "_5"), getattr(row, "_6")\
                            , getattr(row, "_7"), getattr(row, "_8"))
    else:
        newRightBoundary = str(int(getattr(row, "_3")) + 2)
        return nTuple(getattr(row, "_0"), getattr(row, "_1"), "stop_codon"\
                            , getattr(row, "_3"), newRightBoundary\
                            , getattr(row, "_5"), getattr(row, "_6")\
                            , getattr(row, "_7"), getattr(row, "_8"))

# calculate the start codon for a given cds row (ensembl style)
# TODO currently not required, could be exchanged with normal function
def generateStartCodonCDS(row):
    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])
    if getattr(row, "_6") == "+":
        newRightBoundary = str(int(getattr(row, "_3")) + 2)
        return nTuple(getattr(row, "_0"), getattr(row, "_1"), "start_codon"\
                            , getattr(row, "_3"), newRightBoundary\
                            , getattr(row, "_5"), getattr(row, "_6")\
                            , getattr(row, "_7"), getattr(row, "_8"))
    else:
        newLeftBoundary = str(int(getattr(row, "_4") - 2))
        return nTuple(getattr(row, "_0"), getattr(row, "_1"), "start_codon"\
                            , newLeftBoundary, getattr(row, "_4")\
                            , getattr(row, "_5"), getattr(row, "_6")\
                            , getattr(row, "_7"), getattr(row, "_8"))

# calculate the stop codon for a given gene row
def generateStopCodonCDS(row):
    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])
    if getattr(row, "_6") == "+":
        newLeftBoundary = str(int(getattr(row, "_4") + 1))
        newRightBoundary = str(int(getattr(row, "_4") + 3))
        return nTuple(getattr(row, "_0"), getattr(row, "_1"), "stop_codon"\
                            , newLeftBoundary, newRightBoundary\
                            , getattr(row, "_5"), getattr(row, "_6")\
                            , getattr(row, "_7"), getattr(row, "_8"))
    else:
        newLeftBoundary = str(int(getattr(row, "_3")) - 3)
        newRightBoundary = str(int(getattr(row, "_3")) - 1)
        return nTuple(getattr(row, "_0"), getattr(row, "_1"), "stop_codon"\
                            , newLeftBoundary, newRightBoundary\
                            , getattr(row, "_5"), getattr(row, "_6")\
                            , getattr(row, "_7"), getattr(row, "_8"))

def generate_newAnnotation(args, annotationDF, intersectingIDS, features):
    # if codons should be computed, else return the current annotation
    if args.computeCodons:
        newAnnotationDF = pd.DataFrame()
        # calculate start and stop codon with respect to gene, if not existing then take cds instead
        if "gene" in features:
            rows = []
            # using itertuples preserves data types and is usually much faster
            for row in annotationDF.itertuples(index=False, name='Pandas'):
                # only consider IDs existing in the genome file
                if getattr(row, "_0").split(".")[0] not in intersectingIDS:
                    continue
                if getattr(row, "_2").lower() == "gene":
                    rows.append(row)
                    # compute new start and stop codons
                    rows.append(generateStartCodon(row))
                    rows.append(generateStopCodon(row))
                elif getattr(row, "_2").lower() not in ["start_codon", "stop_codon", "end_codon", "begin_codon"]:
                    rows.append(row)

            newAnnotationDF = pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

        else: # cds only
            rows = []
            # using itertuples preserves data types and is usually much faster
            for row in annotationDF.itertuples(index=False, name='Pandas'):
                # only consider IDs existing in the genome file
                if getattr(row, "_0").split(".")[0] not in intersectingIDS:
                    continue
                if getattr(row, "_2").lower() == "cds":
                    rows.append(row)
                    # compute new start and stop codons
                    rows.append(generateStartCodonCDS(row))
                    rows.append(generateStopCodonCDS(row))
                elif getattr(row, "_2").lower() not in ["start_codon", "stop_codon", "end_codon", "begin_codon"]:
                    rows.append(row)

            newAnnotationDF = pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

        return newAnnotationDF
    else:
        # only copy the rows that appear in both annotation and genome files
        rows = []
        for row in annotationDF.itertuples(index=False, name='Pandas'):
            if getattr(row, "_0").split(".")[0] in intersectingIDS:
                rows.append(row)
        newAnnotationDF = pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

        return newAnnotationDF

def validate(args):
    # get the identifiers available in the genome file
    genomeIDs = []
    for line in open(args.genome, "r"):
        if ">" in line:
            genomeIDs.append(re.split('[>. \n]', line)[1])

    # if no identifier is found, something is wrong with the genome file
    if len(genomeIDs) == 0:
        sys.exit("The provided genome file does not contain a valid identifier.\n" \
                +"Ensure that the file is in fasta format and contains atleast one identifier.")

    # Read genome file
    genomeDict = dict()
    with open(args.genome, "r") as genomeFile:
        genomes=genomeFile.read()

    # Split input string at id lines
    splitGenome = re.compile(">.*\n").split(genomes)[1:]
    for i in range(len(genomeIDs)):
        genomeDict[genomeIDs[i]] = splitGenome[i].replace("\n", "")

    # check if the input file is in correct format
    # read the input file
    try:
        annotationDF = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)
    except pd.errors.ParserError as error:
        sys.exit("Error reading the annotation file, please ensure that it is formatted correctly.\n"\
                +"Exiting with:\n   %s" % error)

    # Check number of columns
    numberOfColumns = len(annotationDF.columns)
    if numberOfColumns != 9:
        sys.exit("The annotation file is not conform with the gff standard.\n"\
                +"Expected number of columns: 9, Columns found: %s" % numberOfColumns)

    # get the features and check if gene or cds is available
    features = set(annotationDF[2].unique())
    if len([True for x in features if x.lower() in ["gene", "cds"]]) == 0:
        sys.exit("Invalid annotation file! The annotation file does not contain gene or cds entries.")

    # Check whether the identifiers exist
    genomeIDs = set(genomeIDs)
    annotationIDs = set([a.split(".")[0] for a in annotationDF[0].unique()])

    # handle identifiers:
    # if there is no intersection between annotation file and genome file -> error
    # else: only copy the data for the required identifiers
    intersectingIDs = annotationIDs.intersection(genomeIDs)
    if len(intersectingIDs) == 0:
        sys.exit("Error: The identifiers for genome.fa and annotation.gtf do not match." \
                +" Please, check for empty or missplaced files!")
    else:
        newAnnotationDF = generate_newAnnotation(args, annotationDF, intersectingIDs, features)

    # write new files
    with open(args.annotationOutput, "w") as f:
        f.write("# This annotation file was automatically generated from another annotation file.\n")
    newAnnotationDF.to_csv(args.annotationOutput, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE, mode="a")

    # ids that do only exist in the genome file
    onlyGenome = genomeIDs.difference(annotationIDs)
    # Restructure the .fa file
    newGenomeFile = ""
    for key, value in genomeDict.items():
        if key not in onlyGenome:
            newGenomeFile += ">%s\n" % key
            newGenomeFile += splitString(value, 80)

    # Write output file
    with open(args.genomeOutput, "w") as out:
        out.write(newGenomeFile)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='validate annotation files and generate working files if possible')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True
                        , help="the annotation file to be validated.")
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True
                        , help="the according genome file.")
    parser.add_argument("--annotationOutput", action="store", dest="annotationOutput", required=True
                        , help="output annotation file.")
    parser.add_argument("--genomeOutput", action="store", dest="genomeOutput", required=True
                        , help="output genome file.")
    parser.add_argument("-c", "--computeCodons", action="store_true", dest="computeCodons", default=False
                        , help="flag, if set the start and stop codons for each gene (or alternatively CDS) are calculated."\
                        " If start and stop codons already exist they are recomputed.")
    args = parser.parse_args()

    validate(args)


if __name__ == '__main__':
    main()
