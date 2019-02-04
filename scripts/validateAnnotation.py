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


# This script has 3 major functionalities:
# 1) sanity check the annotation.gtf/gff3 and genome.fa files (default functioniality: required: -a -g)
# 2) create new annotation / genome files containing only intersecting identifiers (--annotationOutput /--genomeOutput)
# 3) (re)compute the start and stop codons for all entries in the annotation file (--annotationOutput -c )
#
# 1) The sanity check     !!! Currently this ignores version numbers of identifiers (as most people seem not to care)!!!
#   - genome.fa      -> check if genome identifiers exist (>XX) if not something is not right
#   - annotation.gtf -> check if the file is readable (by read_csv)
#                    -> check if the file is formatted correctly (9 columns) (might be extendable to check types)
#                    -> check if the file contains "gene" or "cds" features (any typecase)
#
# 2) create new annotation / genome files:
#   - checks which identifiers are common in both files and creates new files containing only those identifiers
#   - allows adding a list of ids which will also be added, despite not being present (--ignore <list> )
#
# 3) (re)compute start / stop codons (requires --annotationOutput to create new file):
#   - if 'gene' feature exists:     compute start/stop codons with respect to the 'gene positions'
#   - if only 'cds' feature exists: compute start/stop codons with respect to the 'cds positions'

sortDict= {"gene": 0, "transcript": 1, "exon": 2, "CDS": 3, "start_codon": 4, "stop_codon": 5}

def cmp_func(a, b):
    if a not in sortDict:
        a = 100
    else:
        a = sortDict[a]

    if b not in sortDict:
        b = 100
    else:
        b = sortDict[b]

    if a < b:
        return -1
    elif a > b:
        return 1
    else:
        return 0

# split the string into strings of given length
def splitString(s, l):
    return "".join([s[i:i+l] + "\n" for i in range(0, len(s), l)])

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
def generateStartCodon(row):
    if getattr(row, "s6") == "+":
        newRightBoundary = str(int(getattr(row, "s3")) + 2)
        return createNTuple(row, s2="start_codon", s4=newRightBoundary)
    else:
        newLeftBoundary = str(int(getattr(row, "s4") - 2))
        return createNTuple(row, s2="start_codon", s3=newLeftBoundary)


# calculate the stop codon for a given cds row (ensembl style)
def generateStopCodon(row):
    if getattr(row, "s2").lower() in ["gene", "exon"]:
        if getattr(row, "s6") == "+":
            newLeftBoundary = str(int(getattr(row, "s4") - 2))
            return createNTuple(row, s2="stop_codon", s3=newLeftBoundary)
        else:
            newRightBoundary = str(int(getattr(row, "s3")) + 2)
            return createNTuple(row, s2="stop_codon", s4=newRightBoundary)
    else:
        if getattr(row, "s6") == "+":
            newLeftBoundary = str(int(getattr(row, "s4") + 1))
            newRightBoundary = str(int(getattr(row, "s4") + 3))
            return createNTuple(row, s2="stop_codon", s3=newLeftBoundary, s4=newRightBoundary)
        else:
            newLeftBoundary = str(int(getattr(row, "s3")) - 3)
            newRightBoundary = str(int(getattr(row, "s3")) - 1)
            return createNTuple(row, s2="stop_codon", s3=newLeftBoundary, s4=newRightBoundary)

# remove all unwanted sequence names from the annotation file and fix the feature notation (typecase)
def cleanse_annotation(annotationDF, intersectingIDs):
    rows = []
    for row in annotationDF.itertuples(index=False, name='Pandas'):
        # skip unwanted sequence names
        if getattr(row, "_0").split(".")[0] not in intersectingIDs:
            continue

        # ensure the right notation of the features
        if getattr(row, "_2").lower() == "cds":
            rows.append(createNTuple(row, s2=getattr(row, "_2").upper()))
        else:
            rows.append(createNTuple(row, s2=getattr(row, "_2").lower()))

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

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

# function to compute missing entries of the annotation (e.g. start / stop codons)
def complete_annotation(cleansedDF):
    # create the geneDict
    geneDict = create_dictionary(cleansedDF)
    print("Completing missing values!")
    # Go through all keys and complete the annotation
    sortedKeys = sorted(geneDict.keys())
    for key in sortedKeys:
        features = set()
        entries = geneDict[key]
        # get the features
        for entry in entries:
            features.add(getattr(entry, "s2"))

        newEntries = []
        for entry in entries:
            newEntries.append(entry)
            # if no CDS create from EXON if possible
            if "CDS" not in features and "exon" in features:
                if getattr(entry, "s2") == "exon":
                    if getattr(entry, "s6") == "+":
                        newEntries.append(createNTuple(entry, s2="CDS", s4=str(int(getattr(entry, "s4") - 3))))
                        features.add("CDS")
                    else:
                        newEntries.append(createNTuple(entry, s2="CDS", s3=str(int(getattr(entry, "s3") + 3))))
                        features.add("CDS")

            # start codon generation
            if "start_codon" not in features:
                if "gene" in features:
                    if getattr(entry, "s2") == "gene":
                        newEntries.append(generateStartCodon(entry))
                        features.add("start_codon")
                elif "exon" in features:
                    if getattr(entry, "s2") == "exon":
                        newEntries.append(generateStartCodon(entry))
                        features.add("start_codon")
                elif "CDS" in features:
                    if getattr(entry, "s2") == "CDS":
                        newEntries.append(generateStartCodon(entry))
                        features.add("start_codon")


            # stop codon generation
            if "stop_codon" not in features:
                if "gene" in features:
                    if getattr(entry, "s2") == "gene":
                        newEntries.append(generateStopCodon(entry))
                        features.add("stop_codon")
                elif "exon" in features:
                    if getattr(entry, "s2") == "exon":
                        newEntries.append(generateStopCodon(entry))
                        features.add("stop_codon")
                elif "CDS" in features:
                    if getattr(entry, "s2") == "CDS":
                        newEntries.append(generateStopCodon(entry))
                        features.add("stop_codon")

        geneDict[key] = newEntries

    return createAnnotationFromDictionary(geneDict)


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
        sys.exit("Invalid annotation file! The annotation file does not contain 'gene' or 'cds' entries.")

    if len([True for x in features if x.lower() in ["exon"]]) == 0:
        sys.exit("No exons detected. Exons can be recomputed using -e (only possible if 'cds' exists)")

    # Check whether the identifiers exist
    genomeIDs = set(genomeIDs)
    annotationIDs = set([a.split(".")[0] for a in annotationDF[0].unique()])

    # ids that do only exist in the genome file / annotation file
    onlyAnnotation = annotationIDs.difference(genomeIDs)
    onlyGenome = genomeIDs.difference(annotationIDs)
    if len(onlyAnnotation) != 0:
        tmpL = ",".join(onlyAnnotation)
        print("Identifiers only present in the annotation file: %s" % tmpL)

    if len(onlyGenome) != 0:
        tmpL = ",".join(onlyGenome)
        print("Identifiers only present in the genome file: %s" % tmpL)

    ignoreA = set([x for x in args.ignoreIDs if x in onlyAnnotation])
    ignoreG = set([x for x in args.ignoreIDs if x in onlyGenome])

    # handle identifiers:
    # if there is no intersection between annotation file and genome file -> error
    # else: only copy the data for the required identifiers
    intersectingIDs = annotationIDs.intersection(genomeIDs)
    if len(intersectingIDs) == 0:
        sys.exit("Error: The identifiers for genome.fa and annotation.gtf do not match." \
                +" Please, check for empty or missplaced files!")

    # if requested, create annotation output
    if args.annotationOutput != "":
        cleansedAnnoDF = cleanse_annotation(annotationDF, intersectingIDs.union(ignoreA))
        newAnnotationDF = complete_annotation(cleansedAnnoDF)

        # if desired sort the dataframe by pos1, pos2, name
        if args.sort:
            print("Sorting annotation file!")
            newAnnotationDF[3] = newAnnotationDF[3].astype(int)
            newAnnotationDF[4] = newAnnotationDF[4].astype(int)
            newAnnotationDF.sort_values([3,4,0], ascending=True, inplace=True)

        # write new files
        with open(args.annotationOutput, "w") as f:
            f.write("# This annotation file was automatically generated from another annotation file.\n")
        newAnnotationDF.to_csv(args.annotationOutput, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE, mode="a")

    # Create newly formatted genome file
    if args.genomeOutput != "":
        # Restructure the .fa file
        newGenomeFile = ""
        for key, value in genomeDict.items():
            if key in intersectingIDs.union(ignoreG):
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
    parser.add_argument("--annotationOutput", action="store", dest="annotationOutput",default="", help="output annotation file.")
    parser.add_argument("--genomeOutput", action="store", dest="genomeOutput", default="", help="output genome file.")
    parser.add_argument("-i", "--ignore", nargs="*", dest="ignoreIDs", default=[], help="IDs that are enforced to be added.")
    parser.add_argument("-s", "--sort", action="store_true", dest="sort", help="sort the file according to 1st position, 2nd position, name")
    args = parser.parse_args()

    validate(args)


if __name__ == '__main__':
    main()
