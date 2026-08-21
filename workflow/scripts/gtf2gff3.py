#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections
import csv
from shutil import copyfile

def generate_dictionaries(args):
    annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    unknown_entries = []

    gene_dict = {}
    cds_dict = {}
    RNA_dict = {}
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        source = getattr(row, "_1")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        score = getattr(row, "_5")
        strand = getattr(row, "_6")
        phase = getattr(row, "_7")
        attributes = getattr(row, "_8")

        if "gene_id " in attributes:
            attribute_list = [x.replace(";", "") for x in list(csv.reader([attributes], delimiter=' ', quotechar='"'))[0]]
            gene_id = attribute_list[attribute_list.index("gene_id")+1]

            if feature.lower() in ["gene", "pseudogene"]:
                gene_dict[gene_id] = (reference_name, source, feature, start, stop, score, strand, phase, attributes)
            elif feature.lower() == "cds":
                cds_dict[gene_id] = (reference_name, source, feature, start, stop, score, strand, phase, attributes)
            elif feature.lower() in ["ncrna", "rrna", "trna", "srna"]:
                RNA_dict[gene_id] = (reference_name, source, feature, start, stop, score, strand, phase, attributes)
            elif feature.lower() in ["transcript"]:
                biotype = ""
                if "gene_biotype" in attribute_list:
                    biotype = attribute_list[attribute_list.index("gene_biotype")+1]
                if biotype in ["ncRNA", "rRNA", "tRNA", "sRNA"]:
                    RNA_dict[gene_id] = (reference_name, source, biotype, start, stop, score, strand, phase, attributes)

            elif feature.lower() not in ["exon", "start_codon", "stop_codon", "transcript"]:
                unknown_entries.append((reference_name, source, feature, start, stop, score, strand, phase, attributes))

        else:
            unknown_entries.append((reference_name, source, feature, start, stop, score, strand, phase, attributes))

    return unknown_entries, gene_dict, cds_dict, RNA_dict

def attributes_to_gff3(attributes, id, locus_tag, parent=""):
    """
    take gff2 format attributes and convert them to gff3
    """

    if "=" in attributes and ";" in attributes:
        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes)]
    else:
        attribute_list = [x.strip(" ") for x in re.split('[;\"]', attributes) if x != ""]

    empty_fields=[]
    for i in range(len(attribute_list)):
        if attribute_list[i] == "":
            empty_fields.extend([i-1, i])
        elif attribute_list[i] == "gene_id":
            empty_fields.extend([i, i+1])
            i+=1

    for index in sorted(empty_fields, reverse=True):
        del attribute_list[index]

    new_attributes="ID=%s;locus_tag=%s;" % (id, locus_tag)
    if parent != "":
        new_attributes += "Parent=%s;" % (parent)

    new_attributes += "".join(["%s=%s;" % (attribute_list[i], attribute_list[i+1]) for i in range(0, len(attribute_list), 2)])

    return new_attributes


def create_gff3_annotation(args):
    unknown_entries, gene_dict, cds_dict, RNA_dict = generate_dictionaries(args)

    union_keys = sorted(list(set().union(gene_dict.keys(), cds_dict.keys(), RNA_dict.keys())))
    nTuple = collections.namedtuple('Pandas', ["seq_name","source","feature","start","stop","score","strand","phase","attribute"])

    cds_count = 1
    gene_count = 1
    RNA_count = 1

    rows = []
    for key in union_keys:
        if key in gene_dict:
            reference_name, source, feature, start, stop, score, strand, phase, attributes = gene_dict[key]

            rows.append(nTuple(reference_name, source, feature, start, stop, score, strand, phase, attributes_to_gff3(attributes, "gene%s" % gene_count, key)))
            if key in cds_dict:
                cds_entries = cds_dict[key]
                rows.append(nTuple(cds_entries[0], cds_entries[1], "CDS", start, stop, cds_entries[5], cds_entries[6], cds_entries[7], attributes_to_gff3(cds_entries[8], "cds%s" % cds_count, key, "gene%s" % gene_count)))
                cds_count += 1

            if key in RNA_dict:
                RNA_entries = RNA_dict[key]
                rows.append(nTuple(RNA_entries[0], RNA_entries[1], RNA_entries[2], start, stop, RNA_entries[5], RNA_entries[6], RNA_entries[7], attributes_to_gff3(RNA_entries[8], "rna%s" % RNA_count, key, "gene%s" % gene_count)))
                RNA_count += 1
            gene_count += 1

        elif key in cds_dict:
            reference_name, source, feature, start, stop, score, strand, phase, attributes = cds_dict[key]

            if strand == "-":
                new_start, new_stop = start - 3, stop
            else:
                new_start, new_stop = start, stop + 3

            rows.append(nTuple(reference_name, source, "gene", new_start, new_stop, score, strand, phase, attributes_to_gff3(attributes, "gene%s"%gene_count, key)))
            rows.append(nTuple(reference_name, source, "CDS", new_start, new_stop, score, strand, phase, attributes_to_gff3(attributes, "cds%s"%cds_count, key, "gene%s"%gene_count)))

            gene_count += 1
            cds_count += 1

        elif key in RNA_dict:
            reference_name, source, feature, start, stop, score, strand, phase, attributes = RNA_dict[key]

            rows.append(nTuple(reference_name, source, "gene", start, stop, score, strand, phase, attributes_to_gff3(attributes, "gene%s"%gene_count, key)))
            rows.append(nTuple(reference_name, source, feature, start, stop, score, strand, phase, attributes_to_gff3(attributes, "rna%s"%rna_count, key, "gene%s"%gene_count)))

            gene_count += 1
            RNA_count += 1

    unknown_count = 1
    for entry in unknown_entries:
        reference_name, source, feature, start, stop, score, strand, phase, attributes = entry
        rows.append(nTuple(reference_name, source, feature, start, stop, score, strand, phase, attributes_to_gff3(attributes, "misc%s" % unknown_count, key)))
        unknown_count += 1

    return pd.DataFrame.from_records(rows, columns=["seqName","source","type","start","stop","score","strand","phase","attribute"])

def which_annotation(args):
    """
    Check if data consitently contains "ID=" or "gene_id ", deciding whether it is in gff2 or gff3 format
    If file has mixed IDs, terminate (for now)
    """
    annotation_df = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    has_gene_id = False
    has_id = False
    for row in annotation_df.itertuples(index=False, name='Pandas'):
        attributes = getattr(row, "_8")

        if "ID=" in attributes:
            has_id = True
        elif "gene_id \"" in attributes:
            has_gene_id = True

        if has_gene_id and has_id:
            print("Mixed gff file, found gene_id and ID identifers. File should be either gff2 or gff3!")
            sys.exit("First occurence of ambiguity: " + row)

    return has_gene_id, has_id

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create gff3 format from gff2 file.')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "input annotation file")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "output annotation file")
    args = parser.parse_args()

    is_gff2, is_gff3 = which_annotation(args)
    if is_gff3:
        copyfile(args.annotation, args.output)
    elif is_gff2:
        annotation_df = create_gff3_annotation(args)
        annotation_df = annotation_df.sort_values(by=["seqName", "start", "stop"])
        with open(args.output, "w") as f:
            f.write("##gff-version 3\n")
        with open(args.output, "a") as f:
            annotation_df.to_csv(f, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
