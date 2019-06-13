#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import csv
import collections

# little helper function to create named tuple without having to always state every argument
def create_ntuple(row, s0=None, s1=None, s2=None, s3=None, s4=None, s5=None, s6=None, s7=None, s8=None):
    nTuple = collections.namedtuple('Pandas', ["s0","s1","s2","s3","s4","s5","s6","s7","s8"])
    try:
        c0, c1, c2 = getattr(row, "_0"), getattr(row, "_1"), getattr(row, "_2")
        c3, c4, c5 = getattr(row, "_3"), getattr(row, "_4"), getattr(row, "_5")
        c6, c7, c8 = getattr(row, "_6"), getattr(row, "_7"), getattr(row, "_8")
    except AttributeError:
        c0, c1, c2 = getattr(row, "s0"), getattr(row, "s1"), getattr(row, "s2")
        c3, c4, c5 = getattr(row, "s3"), getattr(row, "s4"), getattr(row, "s5")
        c6, c7, c8 = getattr(row, "s6"), getattr(row, "s7"), getattr(row, "s8")
    if s0 != None:
        c0 = s0
    if s1 != None:
        c1 = s1
    if s2 != None:
        c2 = s2
    if s3 != None:
        c3 = s3
    if s4 != None:
        c4 = s4
    if s5 != None:
        c5 = s5
    if s6 != None:
        c6 = s6
    if s7 != None:
        c7 = s7
    if s8 != None:
        c8 = s8
    return nTuple(c0, c1, c2, c3, c4, c5, c6, c7, c8)


"""
create a dictionary containing an id for each start,stop,strand combination
"""
def fill_annotation_dict(args):
    # read annotation file
    annotation_df = pd.read_csv(args.annotation, comment="#", sep="\t", header=None)
    annotation_dict = dict()

    # try to guess the file type (gff2/3)
    # number of entries containing ID=
    true_values = sum(annotation_df[8].str.contains("ID=").tolist())
    # number of rows
    nrows = len(annotation_df.index)
    # if 75% of all rows contain ID= it is likely a gff3
    if true_values >= nrows * 0.75:
        # build dictionary
        for row in annotation_df.itertuples(index=False, name='Pandas'):
            # using "gene" for gff3 files as it contains most information
            if getattr(row, "_2").lower() == "gene":
                start = str(getattr(row, "_3"))
                stop = str(getattr(row, "_4"))
                strand = str(getattr(row, "_6"))
                description = getattr(row, "_8")
                attributes = [x.lower() for x in re.split('[;=]', description)]
                # locus_tag
                if "locus_tag" in attributes:
                    locus_tag = attributes[attributes.index("locus_tag") + 1]
                else:
                    locus_tag = attributes[attributes.index("id") + 1]

                # gene name
                if "name" in attributes:
                    name = attributes[attributes.index("name") + 1]
                elif "gene" in attributes:
                    name = attributes[attributes.index("gene") + 1]
                else:
                    name = "NA"

                # update the annotation dictionary
                if start in annotation_dict:
                    if stop in annotation_dict[start]:
                        if strand in annotation_dict[start][stop]:
                            print("Duplicate found!")
                        else:
                            annotation_dict[start][stop][strand] = (locus_tag, name)
                    else:
                        annotation_dict[start][stop] = {}
                        annotation_dict[start][stop][strand] = (locus_tag, name)
                else:
                    annotation_dict[start] = {}
                    annotation_dict[start][stop] = {}
                    annotation_dict[start][stop][strand] = (locus_tag, name)
    else: # gtf
        print("gtf detected")
        # build dictionary
        for row in annotation_df.itertuples(index=False, name='Pandas'):
            # using "cds" for gtf files as it contains most information
            if getattr(row, "_2").lower() == "exon":
                start = str(getattr(row, "_3"))
                stop = str(getattr(row, "_4"))
                strand = str(getattr(row, "_6"))
                description = getattr(row, "_8")
                attributes = [x.lower() for x in re.split('[; ]', description)]
                # locus_tag
                if "locus_tag" in attributes:
                    locus_tag = attributes[attributes.index("locus_tag")+1].replace("\"", "")
                else:
                    locus_tag = attributes[attributes.index("gene_id")+1].replace("\"", "")

                # name
                if "gene_name" in attributes:
                    name = attributes[attributes.index("gene_name")+1].replace("\"", "")
                elif "transcript_name" in attributes:
                    name = attributes[attributes.index("transcript_name")+1].replace("\"", "")
                else:
                    name = "NA"

                # update the annotation dictionary
                if start in annotation_dict:
                    if stop in annotation_dict[start]:
                        if strand in annotation_dict[start][stop]:
                            print("Duplicate found!")
                        else:
                            annotation_dict[start][stop][strand] = (locus_tag, name)
                    else:
                        annotation_dict[start][stop] = {}
                        annotation_dict[start][stop][strand] = (locus_tag, name)
                else:
                    annotation_dict[start] = {}
                    annotation_dict[start][stop] = {}
                    annotation_dict[start][stop][strand] = (locus_tag, name)
    return annotation_dict

"""
Compare the combined gff content to the annotation_dict
"""
def reannotate_ORFs(args):
    annotation_dict = fill_annotation_dict(args)
    # read combined gff
    rows = []
    combined_df = pd.read_csv(args.combinedGFF, comment="#", header=None, sep="\t")
    for row in combined_df.itertuples(index=False, name='Pandas'):
        start = str(getattr(row, "_3"))
        stop = str(getattr(row, "_4"))
        strand = str(getattr(row, "_6"))
        try:
            locus_tag = annotation_dict[start][stop][strand][0]
            name = annotation_dict[start][stop][strand][1]
            rows.append(create_ntuple(row, s8="%s;annotated=%s;name=%s" % (getattr(row, "_8"), locus_tag, name)))
        except KeyError:
            rows.append(row)

    return pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Go over the combinedGFF and annotated ORFs that are already known.')
    parser.add_argument("-c", "--combined", action="store", dest="combinedGFF", required=True, help= "A file in GFF format.")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "The standard annotation file.")
    parser.add_argument("-o", "--output", action="store", dest="outputGFF", required=True, help= "The reannotated output file")
    args = parser.parse_args()

    with open(args.outputGFF, "w") as f:
        f.write("##gff-version 3\n")
    with open(args.outputGFF, "a") as f:
        reannotate_ORFs(args).to_csv(f, header=None, sep="\t", index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
