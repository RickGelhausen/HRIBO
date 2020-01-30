#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import csv
import collections


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
                attributes = getattr(row, "_8")

                if ";" in attributes and "=" in attributes:
                    attribute_list =  [x for x in re.split('[;=]', attributes)]
                else:
                    attribute_list = [x.replace(";", "") for x in list(csv.reader([attributes], delimiter=' ', quotechar='"'))[0]]

                if len(attribute_list) % 2 == 0:
                    for i in range(len(attribute_list)):
                        if i % 2 == 0:
                            attribute_list[i] = attribute_list[i].lower()
                else:
                    print(attributes)
                    sys.exit("Attributes section of gtf/gff is wrongly formatted!")

                # locus_tag
                if "locus_tag" in attribute_list:
                    locus_tag = attribute_list[attribute_list.index("locus_tag") + 1]
                else:
                    locus_tag = attribute_list[attribute_list.index("id") + 1]

                # gene name
                name = ""
                if "name" in attribute_list:
                    name = attribute_list[attribute_list.index("name") + 1]

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
                attributes = [x.replace(";", "") for x in list(csv.reader([description], delimiter=' ', quotechar='"'))[0]]
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
                    name = ""

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
    nTuple = collections.namedtuple('Pandas', ["seq_name","source","feature","start","stop","score","strand","phase","attribute"])
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
            if name == "":
                rows.append(nTuple(getattr(row, "_0"), getattr(row, "_1"), getattr(row, "_2"), start, stop, \
                                   getattr(row, "_5"), strand, getattr(row, "_7"), \
                                   "%s;locus_tag=%s" % (getattr(row, "_8"), locus_tag)))
            else:
                attributes = re.split('[;=]', getattr(row, "_8"))
                attributes[attributes.index("Name")+1] = name
                attributes = ";".join(["%s=%s" % (attributes[i], attributes[i+1]) for i in range(0,len(attributes),2)])
                rows.append(nTuple(getattr(row, "_0"), getattr(row, "_1"), getattr(row, "_2"), start, stop, \
                                   getattr(row, "_5"), strand, getattr(row, "_7"), \
                                   "%s;locus_tag=%s" % (attributes, locus_tag)))
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
        df = reannotate_ORFs(args)
        df.sort_values(by=[0,3,4], inplace=True)
        df.to_csv(f, header=None, sep="\t", index=False, quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
