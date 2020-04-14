#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd
import shlex, subprocess
import collections
import sys
import csv

def call_featureCounts(args):
    """
    set up commandline call for featureCounts, process the featureCounts output
    """
    open(args.output, "w").close()

    bamfiles = sorted(args.bamfiles)
    annotation_df = pd.read_csv(args.annotation, sep="\t", header=None, comment="#")
    features = list(annotation_df[2].unique())

    identifier = "ID" if "ID=" in annotation_df[8][0] else "gene_id"

    tmp_file = os.path.splitext(args.output)[0] + ".tmp"
    commandline_parameters = " -a %s -F GTF -g %s -s %s -T %s -o %s" % (args.annotation, identifier, args.strandness, args.threads, tmp_file)
    if args.assign_to_all:
        commandline_parameters += " -O"
    if args.assign_multi_mappers:
        commandline_parameters += " -M"
    if args.with_fraction:
        commandline_parameters += " --fraction"

    if args.diff_expr:
        header = "," + ",".join([os.path.splitext(os.path.basename(bamfile))[0] for bamfile in bamfiles]) +"\n"
        with open(args.output, "w") as f:
            f.write(header)
        labels = ["s%s" % x for x in range(0, len(bamfiles)+1)]
        nTuple = collections.namedtuple('Pandas', labels)
    else:
        labels = ["s%s" % x for x in range(0, len(bamfiles)+6)]
        nTuple = collections.namedtuple('Pandas', labels)

    for feature in features:
        commandline_call = "featureCounts -t %s " %(feature) + commandline_parameters + " "+ " ".join(bamfiles)
        subprocess_call = shlex.split(commandline_call, posix=False)

        print(commandline_call)
        subprocess.call(subprocess_call)

        try:
            tmp_df = pd.read_csv(tmp_file, skiprows=[1], header=None, sep="\t", comment="#")
            new_rows = []
            for row in tmp_df.itertuples(index=False, name='Pandas'):
                gene_id = getattr(row, "_0")
                chromosome = getattr(row, "_1").split(";")[0]
                start = getattr(row, "_2")
                if type(start) is str:
                    start = start.split(";")[0]

                stop = getattr(row, "_3")
                if type(stop) is str:
                    stop = stop.split(";")[0]
                strand = getattr(row, "_4").split(";")[0]
                length = getattr(row, "_5")
                read_list = [getattr(row, "_%s" % x) for x in range(6, len(row))]

                if args.diff_expr:
                    new_rows.append(nTuple(gene_id, *read_list))
                else:
                    new_rows.append(nTuple(gene_id, chromosome, start, stop, strand, *read_list, feature))

            if args.diff_expr:
                new_df = pd.DataFrame.from_records(new_rows, columns=labels)
                with open(args.output, "a") as f:
                    new_df.to_csv(f, sep=",", header=None, index=False, quoting=csv.QUOTE_NONE)
            else:
                new_df = pd.DataFrame.from_records(new_rows, columns=labels)
                with open(args.output, "a") as f:
                    new_df.to_csv(f, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE)

        except FileNotFoundError:
            sys.exit("temporary file was not found.")

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Call featureCounts and process the output.')
    parser.add_argument("-b", "--bam", nargs="*", dest="bamfiles", required=True, help= "Read sequence files: (.bam)")
    parser.add_argument("-s", "--strandness", action="store", dest="strandness", default=1, help= "Perform strand-specific read counting. 0 (unstranded), 1 (stranded) and 2 (reversely stranded).")
    parser.add_argument("--with_O", action="store_true", dest="assign_to_all", help= "Assign reads to all their overlapping meta-features.")
    parser.add_argument("--with_M", action="store_true", dest="assign_multi_mappers", help= "Multi-mapping reads will also be counted.")
    parser.add_argument("--fraction", action="store_true", dest="with_fraction", help= "Assign fractional counts to features.")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "The annotation to be processed with featureCounts.")
    parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, help= "Number of threads to be used.")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "The output file.")
    parser.add_argument("--for_diff_expr", action="store_true", dest="diff_expr", help= "Number of threads to be used.")
    args = parser.parse_args()

    call_featureCounts(args)

if __name__ == '__main__':
    main()
