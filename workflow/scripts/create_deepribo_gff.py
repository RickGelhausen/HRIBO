#!/usr/bin/env python3
"""Convert DeepRibo predictions to strict GFF3 records."""

import argparse
import csv
import collections

import pandas as pd

import gff_utils


def to_gff3(args):
    inputDF = pd.read_csv(args.predictedORFs, sep=',')
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])

    # extract information from each row and build new dataframe in gff format
    rows = []
    for row in inputDF.itertuples(index=True, name='Pandas'):
        strand = str(getattr(row, "strand"))
        locus = str(getattr(row, "locus"))
        pred = float(getattr(row, "pred"))
        dist = int(getattr(row, "dist"))
        SS_pred_rank = getattr(row, "SS_pred_rank")

        # new content
        chromosome, rest = locus.split(":")
        start, stop = rest.split("-")
        start, stop = int(start), int(stop)

        if SS_pred_rank == 999999:
            continue

        if strand == "+":
            stop += 2
        else:
            start -= 2

        seqName = chromosome
        source = "deepribo"
        feature = "CDS"
        score = pred
        identifier = f"{chromosome}:{start}-{stop}:{strand}"
        attribute = gff_utils.format_attributes(
            [
                ("ID", identifier),
                ("pred_value", str(pred)),
                ("deepribo_distance", str(dist)),
                ("method", "deepribo"),
                ("condition", args.condition),
                ("replicate", args.replicate),
            ]
        )

        # DeepRibo's distance is prediction metadata, not the CDS reading frame.
        # Predictions are complete ORFs beginning at their start codon, so their
        # required CDS phase is zero on either strand.
        rows.append(
            nTuple(
                seqName,
                source,
                feature,
                start,
                stop,
                score,
                strand,
                "0",
                attribute,
            )
        )



    return pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts reperation output to new data frame\
                                     containing specified information and saves it in gff3 format.')
    parser.add_argument("-i", "--inputCSV", action="store", dest="predictedORFs", required=True
                                          , help= "the input file. (created by reparation)")
    parser.add_argument("-c", "--condition", action="store", dest="condition", required=True
                                           , help= "the condition of the current file")
    parser.add_argument("-r", "--replicate", action="store", dest="replicate", required=True
                                           , help= "the condition of the current file")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")

    args = parser.parse_args()

    gff3df = to_gff3(args)
    # Explicit tie-breakers and a stable sort: sorting on score alone left rows
    # of equal score in an order that varied with the pandas build.
    gff3df = gff3df.sort_values(
        by=["score", "seqName", "start", "stop", "strand"], kind="stable"
    )
    with open(args.outputGFF, "w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")
        gff3df.to_csv(
            handle,
            sep="\t",
            header=False,
            index=False,
            quoting=csv.QUOTE_NONE,
        )

if __name__ == '__main__':
    main()
