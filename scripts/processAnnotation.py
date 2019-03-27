import csv
import pandas as pd
import argparse
import re

# little helper function to create named tuple without having to always state every argument
def createNTuple(row, s0=None, s1=None, s2=None, s3=None, s4=None, s5=None, s6=None, s7=None, s8=None):
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


def combineSplit(attrSplit):
    attribute = ""
    attrCount = 0
    openQuotes = False
    for attr in attrSplit:
        if openQuotes:
            if "\"" in attr:
                attribute += attr + "; "
                attrCount = 0
                openQuotes = False
            else:
                attribute += attr + " "

        else:
            if attrCount < 1:
                attrCount += 1
                attribute += attr + " "
            else:
                if attr.count("\"") == 2:
                    attribute += attr + "; "
                    attrCount = 0
                else:
                    attribute += attr + " "
                    openQuotes = True

    return attribute

def process_annotation(args):
    # read input annotation
    annDF = pd.read_csv(args.annotation, sep="\t", comment="#", header=None)

    # check if gff2/gtf or gff3 format
    if annDF[8].str.contains("ID="):
        # gff3
        annDF.to_csv(args.output, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    else:
        # gff2
        # only accept rows that contain gene_id transcript_id combination
        annDF = annDF[annDF[8].str.contains("gene_id") & annDF[8].str.contains("transcript_id")]

        # ensure that the separation of attributes is done correctly
        rows = []
        for row in annDF.itertuples(index=False, name='Pandas'):
            attributes = getattr(row, "_8")
            attrSplit = list(filter(None, re.split('[ ;]', attributes)))

            attributes = combineSplit(attrSplit)
            rows.append(createNTuple(row, s8=attributes))

        annDF = pd.DataFrame.from_records(rows, columns=[0,1,2,3,4,5,6,7,8])

        # write processed annotation to output
        annDF.to_csv(args.output, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='process the annotation for use with plastid')
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "the input annotation file.")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "the processed annotation file.")
    args = parser.parse_args()

    process_annotation(args)


if __name__ == '__main__':
    main()
