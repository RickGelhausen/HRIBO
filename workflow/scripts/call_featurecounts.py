#!/usr/bin/env python
import argparse
import collections
import csv
import os
import shlex
import subprocess
import sys

import pandas as pd


RAW_SCHEMA_VERSION = "#hribo-read-counts-v1"


def bam_labels(bamfiles):
    return [os.path.splitext(os.path.basename(path))[0] for path in bamfiles]


def raw_schema(bamfiles):
    columns = (
        ["Identifier", "Genome", "Start", "Stop", "Strand", "Length"]
        + bam_labels(bamfiles)
        + ["Feature"]
    )
    return RAW_SCHEMA_VERSION + "\t" + "\t".join(columns) + "\n"


def select_features(annotation_df, requested_features, annotation):
    """Return requested feature types that occur in annotation column 3.

    Annotation vocabularies differ between providers.  In particular, many
    bacterial annotations do not use the optional ``sRNA`` type from HRIBO's
    example configuration.  A missing optional type must not abort counting of
    the requested types that are present, but it must remain visible to the
    user so that an unintended spelling or annotation choice is not hidden.
    """
    available_features = list(dict.fromkeys(annotation_df[2].astype(str)))
    if not requested_features:
        return available_features

    features = [
        feature for feature in requested_features if feature in available_features
    ]
    missing_features = [
        feature for feature in requested_features if feature not in available_features
    ]
    if missing_features and features:
        print(
            "WARNING: Skipping requested annotation feature type(s) not present "
            f"in {annotation}: {', '.join(missing_features)}. Feature-type "
            "matching is case-sensitive and uses column 3.",
            file=sys.stderr,
        )

    if not features:
        available = ", ".join(available_features) or "none"
        sys.exit(
            "None of the requested annotation feature types were found in "
            f"{annotation}. Requested: {', '.join(requested_features)}. "
            f"Available column-3 types: {available}. Feature-type matching is "
            "case-sensitive."
        )

    return features


def call_featureCounts(args):
    """
    set up commandline call for featureCounts, process the featureCounts output
    """
    bamfiles = sorted(args.bamfiles, key=lambda s: s.lower())
    try:
        annotation_df = pd.read_csv(
            args.annotation, sep="\t", header=None, comment="#"
        )
    except pd.errors.EmptyDataError:
        annotation_df = pd.DataFrame(columns=range(9))

    if annotation_df.empty:
        if args.features:
            select_features(annotation_df, args.features, args.annotation)
        features = []
    else:
        features = select_features(annotation_df, args.features, args.annotation)

    # Every empty result remains self-describing.  The comment is ignored by
    # legacy pandas readers but tells the mapper how many library columns an
    # empty raw table represents.
    if args.diff_expr:
        header = "Identifier," + ",".join(bam_labels(bamfiles)) + "\n"
    else:
        header = raw_schema(bamfiles)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        handle.write(header)

    if annotation_df.empty:
        return

    # Decided from the whole attribute column rather than its first row: a
    # leading region or source feature without an ID= sent every later lookup to
    # the wrong attribute.
    attribute_column = annotation_df[8].astype(str)
    identifier = "ID" if attribute_column.str.contains("ID=").any() else "gene_id"

    tmp_file = os.path.splitext(args.output)[0] + ".tmp"
    commandline_parameters = [
        "-a",
        args.annotation,
        "-F",
        "GTF",
        "-g",
        identifier,
        "-s",
        str(args.strandness),
        "-T",
        str(args.threads),
        "-o",
        tmp_file,
    ]
    if args.assign_to_all:
        commandline_parameters.append("-O")
    if args.assign_multi_mappers:
        commandline_parameters.append("-M")
    if args.with_fraction:
        commandline_parameters.append("--fraction")

    if args.diff_expr:
        labels = [f"s{x}" for x in range(0, len(bamfiles)+1)]
        nTuple = collections.namedtuple('Pandas', labels)
    else:
        labels = [f"s{x}" for x in range(0, len(bamfiles)+7)]
        nTuple = collections.namedtuple('Pandas', labels)

    for feature in features:
        subprocess_call = [
            "featureCounts",
            "-t",
            feature,
            *commandline_parameters,
            *bamfiles,
        ]

        print(shlex.join(subprocess_call))

        # A stale temporary file from the previous feature would otherwise be
        # read back as though it belonged to this one.
        for leftover in (tmp_file, tmp_file + ".summary"):
            if os.path.exists(leftover):
                os.remove(leftover)

        returncode = subprocess.call(subprocess_call)
        if returncode != 0:
            sys.exit(
                f"featureCounts failed for feature '{feature}' with exit code {returncode}."
            )

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
                read_list = [getattr(row, f"_{x}") for x in range(6, len(row))]

                if args.diff_expr:
                    if feature not in ["gene", "pseudogene"]:
                        new_rows.append(nTuple(f"{chromosome}:{start}-{stop}:{strand}", *read_list))
                else:
                    new_rows.append(nTuple(gene_id, chromosome, start, stop, strand, length, *read_list, feature))

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

    if args.diff_expr:
        df = pd.read_csv(args.output, sep=",")
        df.drop_duplicates(subset=["Identifier"], keep="first", inplace=True)
        with open(args.output, "w") as f:
            df.to_csv(f, sep=",", index=False, quoting=csv.QUOTE_NONE)

    # featureCounts leaves its per-feature table and summary behind.
    for leftover in (tmp_file, tmp_file + ".summary"):
        if os.path.exists(leftover):
            os.remove(leftover)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Call featureCounts and process the output.')
    parser.add_argument("-b", "--bam", nargs="*", dest="bamfiles", required=True, help= "Read sequence files: (.bam)")
    parser.add_argument("-s", "--strandness", action="store", dest="strandness", default=1, help= "Perform strand-specific read counting. 0 (unstranded), 1 (stranded) and 2 (reversely stranded).")
    parser.add_argument("--with_O", action="store_true", dest="assign_to_all", help= "Assign reads to all their overlapping meta-features.")
    parser.add_argument("--with_M", action="store_true", dest="assign_multi_mappers", help= "Multi-mapping reads will also be counted.")
    parser.add_argument("--fraction", action="store_true", dest="with_fraction", help= "Assign fractional counts to features.")
    parser.add_argument("--use_features", nargs="+", dest="features", default=[], help= "The feature type to be used for counting. Default: All features (no gene or pseudgene)")
    parser.add_argument("-a", "--annotation", action="store", dest="annotation", required=True, help= "The annotation to be processed with featureCounts.")
    parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, help= "Number of threads to be used.")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, help= "The output file.")
    parser.add_argument("--for_diff_expr", action="store_true", dest="diff_expr", help= "Number of threads to be used.")
    args = parser.parse_args()

    call_featureCounts(args)


if __name__ == '__main__':
    main()
