#!/usr/bin/env python3
"""Merge duplicate Reparation predictions into a deterministic GFF3 file."""

import argparse
import csv
import collections
import io
import os
import tempfile
from pathlib import Path

import pandas as pd

import gff_utils


GFF3_HEADER = "##gff-version 3\n"


def read_gff(path):
    """Read records while treating zero-byte and header-only GFF3 as empty."""

    try:
        dataframe = pd.read_csv(path, sep="\t", comment="#", header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=range(9))
    if dataframe.shape[1] != 9:
        raise ValueError(f"{path} has {dataframe.shape[1]} columns; expected 9")
    return dataframe


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


def stable_row_key(row):
    return tuple(str(value) for value in row)


def row_probability(row):
    parsed = gff_utils.parse_attributes(getattr(row, "s8"))
    value = gff_utils.first_attribute(parsed, "prob")
    if value == "":
        return float("-inf")
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(
            f"Reparation record {getattr(row, 's0')}:{getattr(row, 's3')}-"
            f"{getattr(row, 's4')} has non-numeric prob={value!r}"
        ) from error


def merged_attributes(identifier, winner, evidence, orf_types):
    """Keep the winning prediction metadata and aggregate shared evidence."""

    winner_pairs = gff_utils.normalize_gff3_attribute_keys(
        gff_utils.split_attributes(getattr(winner, "s8"))
    )
    parsed = gff_utils.parse_attributes(getattr(winner, "s8"))
    name = gff_utils.first_attribute(parsed, "name", default=identifier)
    probability = gff_utils.first_attribute(parsed, "prob")

    replaced = {"id", "name", "orf_type", "evidence", "prob"}
    preserved = [
        (key, value) for key, value in winner_pairs if key.lower() not in replaced
    ]
    pairs = [("ID", identifier), ("Name", name)] + preserved
    pairs.extend(
        [
            ("orf_type", ",".join(sorted(orf_types))),
            ("evidence", " ".join(sorted(evidence))),
            ("prob", probability),
        ]
    )
    return gff_utils.format_attributes(pairs)


def handle_overlap(inputDF):
    """
    read the input gff and merge all duplicate intervals
    """
    nTuple = collections.namedtuple('Pandas', ["seq_name","source","feature","start","stop","score","strand","phase","attribute"])

    # create a dictionary for common ids
    geneDict = create_dictionary(inputDF)

    # run over all entries in the dictionary and combine overlapping ones
    rows = []
    for key in sorted(geneDict):
        # One complete winning record supplies every single-valued scientific
        # field.  Higher probability wins; the full row resolves ties without
        # depending on input order.
        winner = min(
            geneDict[key], key=lambda row: (-row_probability(row), stable_row_key(row))
        )
        evidence = set()
        orftype = set()
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
        attribute = merged_attributes(key, winner, evidence, orftype)

        rows.append(
            nTuple(
                getattr(winner, "s0"),
                "reparation",
                getattr(winner, "s2"),
                getattr(winner, "s3"),
                getattr(winner, "s4"),
                getattr(winner, "s5"),
                getattr(winner, "s6"),
                "0" if str(getattr(winner, "s2")).lower() == "cds" else ".",
                attribute,
            )
        )

    dataframe = pd.DataFrame.from_records(rows, columns=range(9))
    if not dataframe.empty:
        dataframe = dataframe.sort_values(
            by=[0, 3, 4, 6, 2, 8], kind="stable"
        ).reset_index(drop=True)
    return dataframe


def render_gff(dataframe):
    output = io.StringIO()
    output.write(GFF3_HEADER)
    dataframe.to_csv(
        output, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE
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
    parser = argparse.ArgumentParser(description='condense duplicates into one entry')
    parser.add_argument("-i", "--inputGFF", action="store", dest="inputGFF", required=True
                                          , help= "the input file (gff3 format).")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")
    args = parser.parse_args()

    input_dataframe = read_gff(args.inputGFF)
    if input_dataframe.empty:
        newDF = pd.DataFrame(columns=range(9))
    else:
        newDF = handle_overlap(input_dataframe)
    atomic_write(args.outputGFF, render_gff(newDF))


if __name__ == "__main__":
    main()
