#!/usr/bin/env python
"""Summarise the annotation with per-library read counts as a spreadsheet.

One sheet holds every feature and the rest split them by feature type. The
genome lookup, RPKM and translational efficiency arithmetic are shared with the
Reparation and DeepRibo tables in excel_utils.
"""

import argparse

import excel_utils as eu

COLUMNS = (
    eu.identity_columns()
    + eu.LOCUS_COLUMNS
    + eu.MEASURE_COLUMNS
    + eu.SEQUENCE_COLUMNS
    + [("Product", "product"), ("Note", "note")]
)

# Sheet name -> the feature types routed to it, in workbook tab order. The order
# and the names are part of the output people open, so both are fixed here.
# Anything unmatched lands in "miscellaneous", so no feature is silently lost.
SHEETS = [
    ("CDS", ["cds"]),
    ("rRNA", ["rrna"]),
    ("sRNA", ["srna"]),
    ("transcript", ["transcript"]),
    ("5'-UTR", ["5'-utr", "five_prime_utr", "5utr", "five_utr", "5'utr"]),
    ("tRNA", ["trna"]),
    ("pseudogene", ["pseudogene"]),
    ("gene", ["gene"]),
    ("region", ["region"]),
]


def create_excel_file(args):
    context = eu.TableContext(args.genome, args.total_mapped)
    all_df, rows = eu.build_annotation_table(args.reads, context, COLUMNS, source="HRIBO")

    features = [row.feature.lower() for row in rows]

    dataframe_dict = {}
    matched = set()
    for sheet, feature_types in SHEETS:
        matched.update(feature_types)
        selector = [feature in feature_types for feature in features]
        dataframe_dict[sheet] = all_df[selector].reset_index(drop=True)

    dataframe_dict["miscellaneous"] = all_df[
        [feature not in matched for feature in features]
    ].reset_index(drop=True)
    dataframe_dict["all"] = all_df

    eu.excel_writer(args.output_path, dataframe_dict, context.wildcards)


def main():
    parser = argparse.ArgumentParser(
        description="Create an excel file of the annotation with read counts and RPKM."
    )
    parser.add_argument("-g", "--genome", action="store", dest="genome", required=True,
                        help="reference genome")
    parser.add_argument("-t", "--total_mapped_reads", action="store", dest="total_mapped",
                        required=True,
                        help="file containing the total mapped reads for all alignment files.")
    parser.add_argument("-r", "--mapped_reads", action="store", dest="reads", required=True,
                        help="file containing the individual read counts")
    parser.add_argument("-o", "--xlsx", action="store", dest="output_path", required=True,
                        help="output xlsx file")
    args = parser.parse_args()

    create_excel_file(args)


if __name__ == '__main__':
    main()
