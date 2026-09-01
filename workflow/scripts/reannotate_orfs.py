#!/usr/bin/env python
import argparse
import csv
import collections
import io
import os
import tempfile
from pathlib import Path

import pandas as pd

import gff_utils


def generate_annotation_dict(args):
    """
    create dictionary from annotation.
    key : (gene_id, locus_tag, name, gene_name)
    """

    annotation_df = pd.read_csv(args.annotation_path, sep="\t", comment="#", header=None)
    annotation_dict = {}
    gene_dict = {}
    cds_dict = {}

    for row in annotation_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")
        read_list = [getattr(row, "_%s" %x) for x in range(9,len(row))]

        parsed = gff_utils.parse_attributes(attributes)

        new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)

        if feature.lower() == "cds":
            cds_dict[new_key] = (
                gff_utils.first_attribute(parsed, "gene_id", "id"),
                gff_utils.first_attribute(parsed, "locus_tag"),
                gff_utils.first_attribute(parsed, "name", "gene_name"),
                read_list,
                gff_utils.first_attribute(parsed, "old_locus_tag"),
            )
        elif feature.lower() in ["gene", "pseudogene"]:
            gene_dict[new_key] = (
                gff_utils.first_attribute(parsed, "name", "gene_name"),
                gff_utils.first_attribute(parsed, "locus_tag", "gene_id"),
                gff_utils.first_attribute(parsed, "old_locus_tag"),
            )

    for key in cds_dict.keys():
        gene_name = ""
        gene_id, locus_tag, name, read_list, old_locus_tag = cds_dict[key]

        if key in gene_dict:
            gene_name, gene_locus_tag, gene_old_locus_tag = gene_dict[key]

            if locus_tag == "":
                locus_tag = gene_locus_tag
            if old_locus_tag == "":
                old_locus_tag = gene_old_locus_tag

        annotation_dict[key] = (gene_id, locus_tag, name, read_list, gene_name, old_locus_tag)

    return annotation_dict


def reannotate_ORFs(args):
    """
    Compare the combined gff content to the annotation_dict
    """
    nTuple = collections.namedtuple('Pandas', ["seq_name","source","feature","start","stop","score","strand","phase","attribute"])
    annotation_dict = generate_annotation_dict(args)
    # read combined gff
    rows = []
    try:
        combined_df = pd.read_csv(
            args.combinedGFF, comment="#", header=None, sep="\t"
        )
    except pd.errors.EmptyDataError:
        # Zero predictions are a valid result. The caller still writes the GFF3
        # header, giving downstream steps an explicit empty annotation.
        return pd.DataFrame(columns=range(9))
    for row in combined_df.itertuples(index=False, name='Pandas'):
        chromosome = str(getattr(row, "_0"))
        start = int(getattr(row, "_3"))
        stop = int(getattr(row, "_4"))
        strand = str(getattr(row, "_6"))
        key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
        pairs = gff_utils.normalize_gff3_attribute_keys(
            gff_utils.split_attributes(getattr(row, "_8"))
        )
        if key in annotation_dict:
            locus_tag = annotation_dict[key][1]
            name = annotation_dict[key][2]
            gene_name = annotation_dict[key][4]
            old_locus_tag = annotation_dict[key][5]

            # Prefer the gene feature's name, falling back to the feature's own.
            replacement = gene_name or name
            if replacement != "":
                pairs = gff_utils.replace_attribute(pairs, "Name", replacement)
            if locus_tag != "":
                pairs = gff_utils.replace_attribute(pairs, "locus_tag", locus_tag)
            if old_locus_tag != "":
                pairs = gff_utils.replace_attribute(
                    pairs, "old_locus_tag", old_locus_tag
                )

        attributes = gff_utils.format_attributes(pairs)
        feature = getattr(row, "_2")
        phase = "0" if str(feature).lower() == "cds" else getattr(row, "_7")
        rows.append(
            nTuple(
                getattr(row, "_0"),
                getattr(row, "_1"),
                feature,
                start,
                stop,
                getattr(row, "_5"),
                strand,
                phase,
                attributes,
            )
        )

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])


def render_gff(dataframe):
    output = io.StringIO()
    output.write("##gff-version 3\n")
    dataframe.to_csv(
        output, header=None, sep="\t", index=False, quoting=csv.QUOTE_NONE
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
    parser = argparse.ArgumentParser(description='Go over the combinedGFF and annotated ORFs that are already known.')
    parser.add_argument("-c", "--combined", action="store", dest="combinedGFF", required=True, help= "A file in GFF format.")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_path", required=True, help= "The standard annotation file.")
    parser.add_argument("-o", "--output", action="store", dest="outputGFF", required=True, help= "The reannotated output file")
    args = parser.parse_args()

    dataframe = reannotate_ORFs(args)
    dataframe.sort_values(by=[0, 3, 4, 6], inplace=True, kind="stable")
    atomic_write(args.outputGFF, render_gff(dataframe))

if __name__ == '__main__':
    main()
