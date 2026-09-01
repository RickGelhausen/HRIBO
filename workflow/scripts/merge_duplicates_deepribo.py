#!/usr/bin/env python3
"""Merge duplicate DeepRibo predictions into deterministic GFF3 outputs."""

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
GFF_COLUMNS = [
    "seqName",
    "source",
    "type",
    "start",
    "stop",
    "score",
    "strand",
    "phase",
    "attribute",
]
WORK_COLUMNS = GFF_COLUMNS + ["_distance"]


def read_gff(path):
    """Read records while treating zero-byte and header-only GFF3 as empty."""

    try:
        dataframe = pd.read_csv(path, comment="#", sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=range(9))
    if dataframe.shape[1] != 9:
        raise ValueError(f"{path} has {dataframe.shape[1]} columns; expected 9")
    return dataframe


def generate_dictionary(in_df):
    """
    read the input file and create a dictionary containing information on overlapping genes
    (rank, attributes)
    """
    overlap_dict = {}
    for row in in_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        prediction_rank = float(getattr(row, "_5"))
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")
        parsed = gff_utils.parse_attributes(attributes)
        distance = gff_utils.first_attribute(parsed, "deepribo_distance")
        if distance == "":
            # Compatibility with HRIBO intermediates produced before distance
            # metadata moved out of the GFF3 phase column.
            distance = getattr(row, "_7")
        try:
            distance = int(distance)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"DeepRibo record {reference_name}:{start}-{stop}:{strand} "
                "has no integer deepribo_distance"
            ) from error

        key = "%s:%s-%s:%s" % (reference_name, start, stop, strand)
        if key in overlap_dict:
            overlap_dict[key].append((prediction_rank, distance, attributes))
        else:
            overlap_dict[key] = [(prediction_rank, distance, attributes)]

    return overlap_dict


def coordinate_key(key):
    reference_name, interval, strand = key.rsplit(":", 2)
    start, stop = interval.rsplit("-", 1)
    return reference_name, int(start), int(stop), strand

def create_gene_dict(annotation_df):
    gene_dict = {}
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        feature = getattr(row, "_2")
        attributes = getattr(row, "_8")
        if feature == "gene":
            parsed = gff_utils.parse_attributes(attributes)

            gene_dict[gff_utils.first_attribute(parsed, "id")] = (
                gff_utils.first_attribute(parsed, "gene"),
                gff_utils.first_attribute(parsed, "locus_tag"),
                gff_utils.first_attribute(parsed, "old_locus_tag"),
            )
    return gene_dict


def generate_annotation_dict(args):
    annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    parent_dict = create_gene_dict(annotation_df)

    annotation_dict = {}
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")

        if feature not in ["CDS", "cds"]:
            continue

        parsed = gff_utils.parse_attributes(attributes)

        key = "%s:%s-%s:%s" % (reference_name, start, stop, strand)

        parent = gff_utils.first_attribute(parsed, "parent")

        # The feature's own attributes, used where the gene feature has none.
        own_name = gff_utils.first_attribute(parsed, "name", "gene", default=key)
        own_locus_tag = gff_utils.first_attribute(parsed, "locus_tag")
        # Previously read the "locus_tag" attribute here, so old_locus_tag was
        # filled with the current locus tag instead of the old one.
        own_old_locus_tag = gff_utils.first_attribute(parsed, "old_locus_tag")

        if parent in parent_dict:
            name, locus_tag, old_locus_tag = parent_dict[parent]
            if name == "":
                name = own_name
            if locus_tag == "":
                locus_tag = own_locus_tag or "na"
            if old_locus_tag == "":
                old_locus_tag = own_old_locus_tag
        else:
            # old_locus_tag was never initialised on this branch, so it either
            # raised or silently carried over the previous row's value.
            name = own_name
            locus_tag = own_locus_tag
            old_locus_tag = own_old_locus_tag

        annotation_dict[key] = (name, locus_tag, old_locus_tag)

    return annotation_dict

def generate_output_gff(args, overlap_dict):
    """
    write an output file where only the best rank is taken for each overlapping prediction
    """
    nTuple = collections.namedtuple(
        "Pandas",
        [
            "seqName",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
            "distance",
        ],
    )
    annotation_dict = generate_annotation_dict(args)

    rows = []
    rows_plus = []
    for key in sorted(overlap_dict, key=coordinate_key):
        value = overlap_dict[key]
        reference_name, mid, strand = key.rsplit(":", 2)
        start_text, stop_text = mid.rsplit("-", 1)
        start, stop = int(start_text), int(stop_text)
        evidence = set()
        cur_pred_value = -10000.0
        for pred, dist, attribute in value:
            # Handles both the GFF3 and the GTF2 attribute forms.
            parsed = gff_utils.parse_attributes(attribute)

            pred_value = pred
            if pred_value >= cur_pred_value:
                cur_pred_value = pred_value

            # A replicate identifies the evidence precisely; without one, the
            # method and condition are the best that can be said.
            if {"condition", "method", "replicate"} <= parsed.keys():
                evidence.add(parsed["condition"] + "-" + parsed["replicate"])
            elif {"condition", "method"} <= parsed.keys():
                evidence.add(parsed["method"] + "-" + parsed["condition"])

        if key in annotation_dict:
            name, locus_tag, old_locus_tag = annotation_dict[key]
        else:
            name, locus_tag, old_locus_tag = key, "", ""

        new_attribute_pairs = [("ID", key), ("Name", name)]
        if locus_tag != "":
            new_attribute_pairs.append(("locus_tag", locus_tag))
        if old_locus_tag != "":
            new_attribute_pairs.append(("old_locus_tag", old_locus_tag))

        # Distance should agree for duplicate coordinates. If an inconsistent
        # input reaches this point, use the smallest distance deterministically;
        # in particular, preserve the novel-ORF marker (-1) from any replicate.
        distance = min(item[1] for item in value)
        new_attribute_pairs.extend(
            [
                ("pred_value", str(cur_pred_value)),
                ("evidence", " ".join(sorted(evidence))),
                ("deepribo_distance", str(distance)),
            ]
        )
        new_attributes = gff_utils.format_attributes(new_attribute_pairs)

        if cur_pred_value >= 0:
            rows_plus.append(
                nTuple(
                    reference_name,
                    "deepribo",
                    "CDS",
                    start,
                    stop,
                    cur_pred_value,
                    strand,
                    "0",
                    new_attributes,
                    distance,
                )
            )
        rows.append(
            nTuple(
                reference_name,
                "deepribo",
                "CDS",
                start,
                stop,
                cur_pred_value,
                strand,
                "0",
                new_attributes,
                distance,
            )
        )

    return (
        pd.DataFrame.from_records(rows, columns=WORK_COLUMNS),
        pd.DataFrame.from_records(rows_plus, columns=WORK_COLUMNS),
    )


def render_gff(dataframe):
    output = io.StringIO()
    output.write(GFF3_HEADER)
    dataframe.to_csv(
        output, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE
    )
    return output.getvalue()


def stage_output(path, content):
    output = Path(path)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp", text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
        os.chmod(temporary_name, 0o644)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise
    return temporary_name


def atomic_write_outputs(outputs):
    """Replace a related output set together, restoring it on handled errors."""

    outputs = list(outputs)
    targets = [Path(path) for path, _ in outputs]
    if len(set(targets)) != len(targets):
        raise ValueError("output paths must be distinct")
    for target in targets:
        if target.is_dir():
            raise IsADirectoryError(f"output path is a directory: {target}")

    staged = []
    backups = []
    installed = set()
    try:
        for path, content in outputs:
            staged.append((stage_output(path, content), Path(path)))

        # Reserve every backup name before moving any existing output. A failure
        # while staging content or allocating a backup therefore leaves the old
        # pair wholly untouched.
        for target in targets:
            if os.path.lexists(target):
                descriptor, backup = tempfile.mkstemp(
                    dir=target.parent,
                    prefix=f".{target.name}.",
                    suffix=".bak",
                )
                os.close(descriptor)
                os.unlink(backup)
                backups.append((target, backup))
            else:
                backups.append((target, None))

        for target, backup in backups:
            if backup is not None:
                os.replace(target, backup)
        for temporary_name, target in staged:
            os.replace(temporary_name, target)
            installed.add(target)
    except BaseException:
        # Restore in reverse order so a failure replacing the second file cannot
        # leave it paired with the newly installed first file.
        for target, backup in reversed(backups):
            if backup is not None and os.path.lexists(backup):
                os.replace(backup, target)
            elif backup is None and target in installed:
                try:
                    os.unlink(target)
                except FileNotFoundError:
                    pass
        raise
    else:
        for _, backup in backups:
            if backup is not None:
                try:
                    os.unlink(backup)
                except FileNotFoundError:
                    pass
    finally:
        for temporary_name, _ in staged:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def default_plus_output(output_gff):
    output = Path(output_gff)
    if output.suffix.lower() == ".gff":
        return str(output.with_name(f"{output.stem}_plus{output.suffix}"))
    return str(output.with_name(f"{output.name}_plus.gff"))

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='condense duplicates into one entry')
    parser.add_argument("-i", "--inputGFF", action="store", dest="inputGFF", required=True
                                          , help= "the input file (gff3 format).")
    parser.add_argument("-o", "--outputGFF", action="store", dest="outputGFF", required=True
                                           , help= "the output file name (gff3 format)")
    parser.add_argument(
        "-p",
        "--plus-output",
        dest="plus_output",
        help="the non-negative prediction output (default: <output>_plus.gff)",
    )
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True
                                           , help= "annotation file")
    args = parser.parse_args()

    plus_output = args.plus_output or default_plus_output(args.outputGFF)
    if Path(args.outputGFF).resolve() == Path(plus_output).resolve():
        parser.error("--outputGFF and --plus-output must name different files")

    input_dataframe = read_gff(args.inputGFF)
    if input_dataframe.empty:
        newDF = pd.DataFrame(columns=WORK_COLUMNS)
        plusDF = pd.DataFrame(columns=WORK_COLUMNS)
    else:
        orf_dict = generate_dictionary(input_dataframe)
        newDF, plusDF = generate_output_gff(args, orf_dict)
        # Ranks are assigned from this order below, so it has to be reproducible.
        # pandas sorts with an unstable quicksort by default, which ordered rows
        # of equal score differently depending on the pandas and numpy build.
        newDF = newDF.sort_values(
            by=["score", "seqName", "start", "stop", "strand"],
            ascending=[False, True, True, True, True],
            kind="stable",
        )
        newDF = newDF.reset_index()
        counter = 1
        updated_attributes = []
        for distance, attributes in zip(newDF["_distance"], newDF["attribute"]):
            # Some entries are neither 0 nor -1. They retain the historical
            # non-novel sentinel because DeepRibo does not define them as novel.
            if distance == -1:
                novel_rank = counter
                counter += 1
            else:
                novel_rank = 999999
            pairs = gff_utils.split_attributes(attributes)
            pairs.append(("novel_rank", str(novel_rank)))
            updated_attributes.append(gff_utils.format_attributes(pairs))

        newDF["attribute"] = updated_attributes

        newDF["score"] = newDF.index + 1
        newDF = newDF.drop(columns=["index"])

        plusDF = plusDF.sort_values(
            by=["seqName", "start", "stop", "strand"], kind="stable"
        )

    newDF = newDF.drop(columns=["_distance"])
    plusDF = plusDF.drop(columns=["_distance"])

    atomic_write_outputs(
        [
            (args.outputGFF, render_gff(newDF)),
            (plus_output, render_gff(plusDF)),
        ]
    )

if __name__ == "__main__":
    main()
