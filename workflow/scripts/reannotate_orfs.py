#!/usr/bin/env python
import argparse
import os
import pandas as pd
import csv
import collections
import sys

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
    combined_df = pd.read_csv(args.combinedGFF, comment="#", header=None, sep="\t")
    for row in combined_df.itertuples(index=False, name='Pandas'):
        chromosome = str(getattr(row, "_0"))
        start = str(getattr(row, "_3"))
        stop = str(getattr(row, "_4"))
        strand = str(getattr(row, "_6"))
        key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
        try:
            locus_tag = annotation_dict[key][1]
            name = annotation_dict[key][2]
            gene_name = annotation_dict[key][4]
            old_locus_tag = annotation_dict[key][5]

            pairs = gff_utils.split_attributes(getattr(row, "_8"))

            # Prefer the gene feature's name, falling back to the feature's own.
            replacement = gene_name or name
            if replacement != "":
                pairs = gff_utils.replace_attribute(pairs, "Name", replacement)
            attributes = gff_utils.format_attributes(pairs)

            if locus_tag != "":
                attributes += "locus_tag=%s;" % locus_tag
            if old_locus_tag != "":
                attributes += "old_locus_tag=%s;" % old_locus_tag

            rows.append(nTuple(getattr(row, "_0"), getattr(row, "_1"), getattr(row, "_2"), start, stop, \
                               getattr(row, "_5"), strand, getattr(row, "_7"), attributes))

        except KeyError:
            rows.append(row)

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Go over the combinedGFF and annotated ORFs that are already known.')
    parser.add_argument("-c", "--combined", action="store", dest="combinedGFF", required=True, help= "A file in GFF format.")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation_path", required=True, help= "The standard annotation file.")
    parser.add_argument("-o", "--output", action="store", dest="outputGFF", required=True, help= "The reannotated output file")
    args = parser.parse_args()

    with open(args.outputGFF, "w") as f:
        f.write("##gff-version 3\n")
    with open(args.outputGFF, "a") as f:
        df = reannotate_ORFs(args)
        df.sort_values(by=[0, 3, 4, 6], inplace=True, kind="stable")
        df.to_csv(f, header=None, sep="\t", index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
