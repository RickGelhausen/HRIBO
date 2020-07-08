#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import csv
import collections
import sys


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

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attribute_list)
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")

        if feature.lower() == "cds":
            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            name = ""
            if "name" in attribute_list:
                name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                name = attribute_list[attribute_list.index("gene_name")+1]

            gene_id = ""
            if "gene_id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("gene_id")+1]
            elif "id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("id")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            cds_dict[new_key] = (gene_id, locus_tag, name, read_list, old_locus_tag)
        elif feature.lower() in ["gene","pseudogene"]:
            gene_name = ""
            if "name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("gene_name")+1]

            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]
            elif "gene_id" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("gene_id")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            gene_dict[new_key] = (gene_name, locus_tag, old_locus_tag)

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

            attribute_list = [x.strip(" ") for x in re.split('[;=]', getattr(row, "_8")) if x != ""]

            if gene_name != "":
                attribute_list[attribute_list.index("Name")+1] = gene_name
                attributes = ";".join(["%s=%s" % (attribute_list[i], attribute_list[i+1]) for i in range(0,len(attribute_list),2)]) +";"
            elif name != "":
                attribute_list[attribute_list.index("Name")+1] = name
                attributes = ";".join(["%s=%s" % (attribute_list[i], attribute_list[i+1]) for i in range(0,len(attribute_list),2)]) +";"
            else:
                attributes = ";".join(["%s=%s" % (attribute_list[i], attribute_list[i+1]) for i in range(0,len(attribute_list),2)]) +";"

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
        df.sort_values(by=[0,3,4], inplace=True)
        df.to_csv(f, header=None, sep="\t", index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
