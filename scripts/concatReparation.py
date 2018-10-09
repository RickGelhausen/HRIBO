#!/usr/bin/env python
'''This script takes input files generated by
reparationGFF and merges all replicates into one file, for each condition.
'''

import os
import argparse
import pandas as pd

def concatGff(gffFiles, prefix, inputFolder):
    # create dataframes of all non-empty files
    dataFrames = []
    for file in gffFiles:
        if os.stat(os.path.join(inputFolder,file)).st_size != 0:
            dataFrames.append(pd.read_csv(os.path.join(inputFolder,file), sep='\t', skiprows=[0], header=None))

    # check if dataframe exist for concatination
    if len(dataFrames) != 0:
        mergedGff = pd.concat(dataFrames)
        outputFile = os.path.join(inputFolder,"%s.reparation.gff" %prefix)
        ### Handling output
        # Setting up result file with header
        with open(outputFile, 'w') as f:
            f.write("##gff-version 3\n")

        # write to file
        with open(outputFile, 'a') as f:
            mergedGff.to_csv(f, sep="\t", header=False, index=False)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts reperation output to new data frame\
                                     containing specified information and saves it in gff3 format.')
    parser.add_argument("-i", "--inputFolder", action="store", dest="tracks", required=True
                                             , help= "the folder containing the output of reparation.")
    parser.add_argument("-d", "--delete", action="store_true", dest="delete", default=False
                                        , help="delete the initial files after merging is complete")
    args = parser.parse_args()

    # get all file names
    gffFiles = []
    for file in os.listdir(args.tracks):
        if file.endswith(".reparation.gff"):
            gffFiles.append(file)

    # create list of prefixes
    prefixList = set([x.split("-")[0] for x in gffFiles])

    # Merge all files with the same prefix
    for prefix in prefixList:
        concatGff(sorted([x for x in gffFiles if ("-" in x) and (x.split("-")[0] == prefix)]), prefix, args.tracks)

    # delete the initial files
    if args.delete:
        for file in gffFiles:
            os.remove(os.path.join(args.tracks,file))

if __name__ == '__main__':
    main()