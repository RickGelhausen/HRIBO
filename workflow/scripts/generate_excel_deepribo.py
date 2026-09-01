#!/usr/bin/env python
"""Summarise the DeepRibo ORF predictions as a spreadsheet.

The table layout is the only thing specific to this script; the genome lookup,
RPKM and translational efficiency arithmetic are shared with the Reparation and
annotation tables in excel_utils.
"""

import argparse

import excel_utils as eu

COLUMNS = (
    eu.identity_columns()
    + [
        # Named as in the overview table, so the same quantity has one name
        # across the workbooks. DeepRibo puts its prediction rank in the score
        # column. Novel rank is prediction metadata stored in a lowercase GFF3
        # attribute; excel_utils still reads legacy phase-based files.
        ("Deepribo_score", "pred_value"),
        ("Deepribo_rank", "score"),
        ("Novel_rank", "novel_rank"),
    ]
    + eu.LOCUS_COLUMNS
    + eu.MEASURE_COLUMNS
    + [("Evidence", "evidence")]
    + eu.SEQUENCE_COLUMNS
)


def create_excel_file(args):
    context = eu.TableContext(args.genome, args.total_mapped)
    cds_df, _ = eu.build_annotation_table(args.reads, context, COLUMNS)
    cds_df = cds_df.sort_values(by=["Genome", "Start", "Stop", "Strand"], kind="stable")

    eu.excel_writer(args.output_path, {"CDS": cds_df}, context.wildcards)


def main():
    parser = argparse.ArgumentParser(
        description="Create an excel file of the DeepRibo predictions with read counts and RPKM."
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
