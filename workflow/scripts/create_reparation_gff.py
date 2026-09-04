#!/usr/bin/env python3
"""Convert Reparation predictions to strict GFF3 records."""

import argparse
import csv
import collections

import pandas as pd

import gff_utils


def createNTuple(args, row):
    nTuple = collections.namedtuple('Pandas', ["seqName","source","type","start","stop","score","strand","phase","attribute"])
    # txt file content
    ORF_locus = getattr(row, "ORF_locus")
    strand = getattr(row, "strand")
    length = str(getattr(row, "length"))
    ribo_count = str(getattr(row, "ribo_count"))
    ribo_rpkm = str(getattr(row, "ribo_rpkm"))
    ribo_coverage = str(getattr(row, "ribo_coverage"))
    SD_score = str(getattr(row, "SD_score"))
    SD_pos = str(getattr(row, "SD_pos"))
    prob = str(getattr(row, "prob"))
    ORF_type = getattr(row, "ORF_type")
    Reference = str(getattr(row, "Reference"))
    Distance_from_aTIS = str(getattr(row, "Distance_from_aTIS"))

    # new content
    # Sequence identifiers may themselves contain colons.  REPARATION appends
    # the numeric coordinate suffix after the final colon.
    chromosome, rest = ORF_locus.rsplit(":", 1)
    start, stop = rest.split("-")
    # modify coordinates to include stop codon
    if strand == '+':
       start = str(int(start)) # due to bug in reparation
       stop = str(int(stop) + 3)
    if strand == '-':
       start = str(int(start) - 3)
       stop = str(int(stop)) # due to bug in reparation

    seqName = chromosome
    source = "reparation"
    type = "CDS"
    score = "."
    # Reparation emits complete ORFs beginning at a start codon. GFF3 therefore
    # requires CDS phase zero on both strands.
    phase = "0"
    identifier = f"{chromosome}:{start}-{stop}:{strand}"
    attribute = gff_utils.format_attributes(
        [
            ("ID", identifier),
            ("Name", identifier),
            ("orf_type", ORF_type),
            ("length", length),
            ("ribo_count", ribo_count),
            ("ribo_rpkm", ribo_rpkm),
            ("ribo_coverage", ribo_coverage),
            ("sd_score", SD_score),
            ("sd_pos", SD_pos),
            ("prob", prob),
            ("reference", Reference),
            ("distance_from_atis", Distance_from_aTIS),
            ("condition", args.condition),
            ("replicate", args.replicate),
            ("method", "reparation"),
        ]
    )

    return nTuple(seqName, source, type, start, stop, score, strand, phase, attribute)


def to_gff3(args):
    inputDF = pd.read_csv(args.predictedORFs, sep='\t')

    # extract information from each row and build new dataframe in gff format
    rows = []
    for row in inputDF.itertuples(index=True, name='Pandas'):
        rows.append(createNTuple(args, row))

    return pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts reperation output to new data frame\
                                     containing specified information and saves it in gff3 format.')
    parser.add_argument("-i", "--inputTXT", action="store", dest="predictedORFs", required=True
                                          , help= "the input file. (created by reparation)")
    parser.add_argument("-c", "--condition", action="store", dest="condition", required=True
                                          , help= "the condition of the current file")
    parser.add_argument("-r", "--replicate", action="store", dest="replicate", required=True
                                          , help= "the replicate of the current file")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")

    args = parser.parse_args()
    gff3df = to_gff3(args)

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
