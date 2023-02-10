#!/usr/bin/env python
import argparse
import pandas as pd
import collections
import excel_utils as eu

def get_unique(in_list):
    seen = set()
    seen_add = seen.add
    return [x for x in in_list if not (x in seen or seen_add(x))]

def parse_orfs(args):

    with open(args.total_mapped, "r") as f:
        total = f.readlines()

    wildcards = []
    for line in total:
        wildcard, reference_name, value = line.strip().split("\t")
        wildcards.append(wildcard)

    wildcards = get_unique(wildcards)
    #read bed file
    read_df = pd.read_csv(args.reads, comment="#", header=None, sep="\t")

    # read gff file
    main_sheet = []

    header = ["Orientation","Class", "Feature count"] + wildcards
    prefix_columns = len(read_df.columns) - len(wildcards)
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    decode = {"ncrna" : "ncRNA", "srna" : "sRNA", "5'-utr" : "5'-UTR", "cds" : "CDS", "rrna" : "rRNA", "trna" : "tRNA", "transcript" : "transcript", "pseudogene" : "pseudogene", "total" : "total"}
    feature_list = ["ncrna", "srna", "5'-utr", "cds", "rrna", "trna", "transcript", "pseudogene", "total"]

    alias_dict = {"5'utr" : "5'-utr", "five_prime_utr" : "5'-utr", "5utr" : "5'-utr", "3'utr" : "3'-utr", "three_prime_utr" : "3'-utr", "3utr" : "3'-utr"}
    read_dict = collections.OrderedDict()
    count_dict = collections.OrderedDict()

    for f in feature_list:
        read_dict[f] = [0] * len(wildcards)
        count_dict[f] = 0

    main_sheet = []
    for row in read_df.itertuples(index=False, name='Pandas'):
        reference_name = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        feature = getattr(row, "_2")

        read_list = [getattr(row, "_%s" %x) for x in range(prefix_columns,len(row))]
        if feature.lower() in read_dict:
            for idx, value in enumerate(read_list):
                read_dict[feature.lower()][idx] += value
                read_dict["total"][idx] += value
            count_dict[feature.lower()] += 1
            count_dict["total"] += 1
        elif feature.lower() in alias_dict:
            if alias_dict[feature.lower()] in read_dict:
                for idx, value in enumerate(read_list):
                    read_dict[alias_dict[feature.lower()]][idx] += value
                    read_dict["total"][idx] += value
                count_dict[alias_dict[feature.lower()]] += 1
                count_dict["total"] += 1
        else:
            print("feature not usable: " + feature)

    for key, val in read_dict.items():
        result = ["sense", decode[key], count_dict[key]] + [float("%.2f" % v) for v in val]

        main_sheet.append(nTuple(*result))

    main_df = pd.DataFrame.from_records(main_sheet, columns=[header[x] for x in range(len(header))])

    dataframe_dict = { "Main" : main_df }

    eu.excel_writer(args.output, dataframe_dict, wildcards)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create excel table with read information.')
    parser.add_argument("-r", "--mapped_reads", action="store", dest="reads", required=True, help= "file containing the individual read counts")
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped", required=True, help= "total mapped reads file")
    parser.add_argument("-o", "--xlsx", action="store", dest="output", required=True, help= "output xlsx file")
    args = parser.parse_args()

    parse_orfs(args)

if __name__ == '__main__':
    main()
