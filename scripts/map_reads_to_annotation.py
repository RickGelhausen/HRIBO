#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections
import csv


def read_file_to_dictionary(args):
    """
    read feature counts read file into dictionary
    """
    read_df = pd.read_csv(args.input, sep="\t", comment="#", header=None)

    read_dict = {}
    for row in read_df.itertuples(index=False, name='Pandas'):
        gene_id = getattr(row, "_0")
        reference_name = getattr(row, "_1")
        start = getattr(row, "_2")
        stop = getattr(row, "_3")
        strand = getattr(row, "_4")
        feature = getattr(row, "_%s" % (len(read_df.columns)-1))

        if ";" in str(start):
            reference_name = reference_name.split(";")
            start = start.split(";")
            stop = stop.split(";")
            strand = strand.split(";")

            for idx in range(len(start)):
                key = (feature, reference_name[idx], start[idx], stop[idx], strand[idx])
                value = []
                for idx in range(6, len(read_df.columns)-1):
                    value.append(getattr(row, "_%s" % idx))

                read_dict[key] = value

        # _6 +
        else:
            key = (feature, reference_name, start, stop, strand)
            value = []
            for idx in range(6, len(read_df.columns)-1):
                value.append(getattr(row, "_%s" % idx))

            read_dict[key] = value

    return read_dict, len(read_df.columns) - 7


def map_reads_to_annotation(args):
    """
    map the reads in the dictionary to the annotation
    """

    annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)
    read_dict, read_number = read_file_to_dictionary(args)
    read_size = read_number
    read_number += len(annotation_df.columns)
    name_list = ["s%s" % str(x) for x in range(read_number)]
    nTuple = collections.namedtuple('Pandas', name_list)

    rows = []
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        info = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        score = getattr(row, "_5")
        strand = getattr(row, "_6")
        phase = getattr(row, "_7")
        attributes = getattr(row, "_8")

        key = (feature, reference_name, start, stop, strand)
        try:
            result = [reference_name, info, feature, start, stop, score, strand, phase, attributes] + read_dict[key]
        except KeyError:
            result = [reference_name, info, feature, start, stop, score, strand, phase, attributes] + [0]*read_size

        rows.append(nTuple(*result))

    return pd.DataFrame.from_records(rows, columns=[x for x in range(len(name_list))])

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='map reads from feature counts to the original annotation')

    parser.add_argument("-i", "--input", action="store", dest="input", required=True, help= "raw featurecounts formats.")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "annotation the reads will be mapped to.")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "output gtf")
    args = parser.parse_args()

    map_reads_to_annotation(args).to_csv(args.output, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)




if __name__ == '__main__':
    main()
