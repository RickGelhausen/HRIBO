#!/usr/bin/env python
import argparse
import collections
import csv
import io
import os
import tempfile
from pathlib import Path

import pandas as pd


RAW_SCHEMA_VERSION = "#hribo-read-counts-v1"
MAPPED_SCHEMA_VERSION = "#hribo-gff-read-counts-v1"
GFF_COLUMNS = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]


def read_schema(path):
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(RAW_SCHEMA_VERSION + "\t"):
                return line.rstrip("\r\n").split("\t")[1:]
    return []


def read_file_to_dictionary(args):
    """
    read feature counts read file into dictionary
    """
    schema = read_schema(args.input)
    try:
        read_df = pd.read_csv(args.input, sep="\t", comment="#", header=None)
    except pd.errors.EmptyDataError:
        read_count = max(len(schema) - 7, 0)
        library_names = schema[6:-1] if schema else []
        return {}, read_count, library_names

    read_dict = {}
    for row in read_df.itertuples(index=False, name='Pandas'):
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

    read_count = len(read_df.columns) - 7
    library_names = schema[6:-1] if schema else []
    if len(library_names) != read_count:
        library_names = [f"read_count_{index + 1}" for index in range(read_count)]
    return read_dict, read_count, library_names


def map_reads_to_annotation(args):
    """
    map the reads in the dictionary to the annotation
    """

    try:
        annotation_df = pd.read_csv(
            args.annotation, sep="\t", comment="#", header=None
        )
    except pd.errors.EmptyDataError:
        annotation_df = pd.DataFrame(columns=range(9))
    read_dict, read_number, library_names = read_file_to_dictionary(args)
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

    return (
        pd.DataFrame.from_records(rows, columns=[x for x in range(len(name_list))]),
        library_names,
    )


def render_mapped_reads(dataframe, library_names):
    output = io.StringIO()
    # Mapped read tables extend the nine GFF columns with one count per library,
    # so they are deliberately described by an HRIBO schema comment rather than
    # falsely claiming to be strict nine-column GFF3.
    output.write(
        MAPPED_SCHEMA_VERSION
        + "\t"
        + "\t".join(GFF_COLUMNS + library_names)
        + "\n"
    )
    dataframe.to_csv(
        output, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE
    )
    return output.getvalue()


def atomic_write(path, content):
    output = Path(path)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp", text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
        os.chmod(temporary_name, 0o644)
        os.replace(temporary_name, output)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='map reads from feature counts to the original annotation')

    parser.add_argument("-i", "--input", action="store", dest="input", required=True, help= "raw featurecounts formats.")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "annotation the reads will be mapped to.")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "output gtf")
    args = parser.parse_args()

    dataframe, library_names = map_reads_to_annotation(args)
    atomic_write(args.output, render_mapped_reads(dataframe, library_names))




if __name__ == '__main__':
    main()
