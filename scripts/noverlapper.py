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


def getLongestInterval(intervals):
    maxi = -1
    # sort intervals by interval ends
    intervals.sort(key=itemgetter(1))
    # create list of weights
    weights = [interval[2] for interval in intervals]

    # find longest non-overlapping interval
    for i in range(1,len(intervals)):
        for j in range(0,i):
            # if non-overlapping: calculate max weights
            if intervals[j][1] <= intervals[i][0]:
                weights[i]=max(weights[i],weights[j]+intervals[i][2])
            maxi = max(weights[i], maxi)

    # Traceback
    optimalInterval=[]
    for i in range(len(intervals)-1,-1,-1):
        if weights[i] == maxi:
            maxi -= intervals[i][2]
            optimalInterval.append(intervals[i])

    return optimalInterval[::-1]

def handleOverlap(args):
    # read input file into dataframe
    inputDF = pd.read_csv(args.inputGFF, sep='\t', header=None)

    # remove entries with duplicate ranges
    inputDF.drop_duplicates(subset=[0,3,4], keep="first",inplace=True)

    # get different identifiers
    identifiers = set(inputDF[0])

    finalTuples = []
    # for each identifier find the longest non-overlapping interval
    for id in identifiers:
        # create a list of tuples used in the following algorithm
        rowTuples = [(row[3], row[4], row[4]-row[3], i) for i, row in inputDF.loc[inputDF[0] == id].iterrows()]

        # find the longest non-overlapping intervalc
        finalTuples += getLongestInterval(rowTuples)

    # create new dataframe containing only the non-overlapping entries for each identifier
    finalTuples.sort(key=itemgetter(3))
    newDF = inputDF.loc[[x[3] for x in finalTuples],:]
    with open(args.outputGFF, 'w') as f:
        newDF.to_csv(f, sep="\t", header=False, index=False)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='removes duplicate intervals and finds\
                                                    the longest non-overlapping interval')
    parser.add_argument("-i", "--inputGFF", action="store", dest="inputGFF", required=True
                                          , help= "the input file (gff3 format).")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")
    args = parser.parse_args()

    handleOverlap(args)


if __name__ == '__main__':
    main()
