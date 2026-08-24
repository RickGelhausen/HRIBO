#!/usr/bin/env python
'''This script takes input gff3 files and handles
overlapping intervals, by merging duplicates.
'''
import pandas as pd
import argparse
import os
import csv
import collections

import gff_utils

def create_dictionary(inputDF):
    """
    function to create a dictionary {geneID : namedtuple}
    """
    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])

    geneDict = dict()
    for row in inputDF.itertuples(index=False, name='Pandas'):
        parsed = gff_utils.parse_attributes(getattr(row, "_8"))
        if "id" in parsed:
            geneID = parsed["id"]
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
        cur_highest_proba = 0
        for row in geneDict[key]:
            parsed = gff_utils.parse_attributes(getattr(row, "s8"))

            # A replicate identifies the evidence precisely; without one, the
            # method and condition are the best that can be said.
            if {"condition", "method", "replicate"} <= parsed.keys():
                evidence.add(
                    parsed["method"] + "-" + parsed["condition"] + "-" + parsed["replicate"]
                )
            elif {"condition", "method"} <= parsed.keys():
                evidence.add(parsed["method"] + "-" + parsed["condition"])

            if "orf_type" in parsed:
                orftype.add(parsed["orf_type"])
            if "prob" in parsed:
                cur_proba = float(parsed["prob"])

            if cur_proba > cur_highest_proba:
                cur_highest_proba = cur_proba

        # sorted(): joining a set directly makes the output depend on the process
        # hash seed, so the same input produced a different file on every run.
        attribute = "ID="+key+";Name="+key+";ORF_type="+",".join(sorted(orftype))+";Evidence="+" ".join(sorted(evidence))+";Prob=" + str(cur_highest_proba)

        rows.append(nTuple(getattr(sampleRow, "s0"),"reparation", getattr(sampleRow, "s2"), getattr(sampleRow, "s3"), \
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
