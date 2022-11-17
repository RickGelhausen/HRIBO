#!/usr/bin/env python

import argparse
from pathlib import Path

import lib.io as io
import lib.misc as misc
import lib.annotation as ann
import lib.metagene as mg
import lib.plotting as plotting

from lib.alignment import IntervalReader

def create_metagene_figures(start_coverage_dict, stop_coverage_dict, read_length_list, meta_dir, mapping_method, normalization_method, positions_out_ORF, positions_in_ORF, color_list):
    """
    Create metagene profiles for all chromosomes and strands for a given mapping and normalization method.
    """

    start_coverage_dict, stop_coverage_dict = misc.equalize_dictionary_keys(start_coverage_dict, stop_coverage_dict, positions_out_ORF, positions_in_ORF)

    df_start_dict = misc.create_data_frame(start_coverage_dict, positions_out_ORF, positions_in_ORF, "start")
    df_stop_dict = misc.create_data_frame(stop_coverage_dict, positions_out_ORF, positions_in_ORF, "stop")

    window_size = positions_out_ORF + positions_in_ORF

    fig_list = []
    for chromosome in df_start_dict:
        df_start = df_start_dict[chromosome]
        df_stop = df_stop_dict[chromosome]

        if normalization_method == "window":
            df_start = misc.window_normalize_df(df_start, window_size)
            df_stop = misc.window_normalize_df(df_stop, window_size)

        max_y = None
        tmp_fig, max_y = plotting.plot_metagene_profiles(df_start, df_stop, read_length_list, f"<b>{chromosome}</b>", max_y=max_y, color_list=color_list)

        fig_list.append((chromosome, mapping_method, tmp_fig))

    io.create_excel_file(df_start_dict, meta_dir / f"{mapping_method}_readcounts_start.xlsx")
    io.create_excel_file(df_stop_dict, meta_dir / f"{mapping_method}_readcounts_stop.xlsx")

    return fig_list

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Perform metagene profiling analysis.", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-b", "--alignment_file_path", action="store", dest="alignment_file_path", type=Path, required=True\
                                                    , help="The path to the alignment (.sam or .bam) file\
                                                            preferably in the format <METHOD>-<CONDITION>-<REPLICATE>(.sam|.bam).")
    parser.add_argument("-a", "--annotation_file_path", action="store", dest="annotation_file_path", type=Path, required=True\
                                                      , help="The path to the output directory.")
    parser.add_argument("-g", "--genome_file_path", action="store", dest="genome_file_path", type=Path, required=True\
                                                  , help="The path to the genome file.")
    parser.add_argument("-o", "--output_dir_path", action="store", dest="output_dir_path", type=Path, required=True\
                                                 , help="The path to the output directory.")
    parser.add_argument("-r", "--read_lengths", action="store", dest="read_lengths", type=str, default="25-34"\
                                              , help="The read lengths to be considered for the metagene-profiling.")
    parser.add_argument("-m", "--mapping_methods", nargs="+", action="store", dest="mapping_methods", default=["fiveprime,threeprime"]\
                                                 , help="The mapping method used for the alignment.")
    parser.add_argument("-n", "--normalization_methods", nargs="+", action="store", dest="normalization_methods", default=["raw"]\
                                                       , help="The normalization method used for the read counts (raw, cpm, window). Default: raw.")
    parser.add_argument("--length_cutoff", action="store", dest="length_cutoff", type=int, default=50\
                                         , help="The minimum length of an ORF to be considered for the metagene-profiling.")
    parser.add_argument("--neighboring_genes_distance", action="store", dest="neighboring_genes_distance", type=int, default=50\
                                            , help="The distance to check for overlapping genes. Default: 50.")
    parser.add_argument("--filtering_methods", nargs="+", dest="filtering_methods", default=["overlap", "rpkm", "length"]\
                                             , help="The filtering methods to be used for the annotation filtering. Default: overlap, rpkm, length.")
    parser.add_argument("--rpkm_threshold", action="store", dest="rpkm_threshold", type=float, default=10.0\
                                          , help="The RPKM threshold to filter genes. Default: 10.0.")
    parser.add_argument("--color_list", nargs="+", dest="color_list", required=False\
                                      , default=[], help="List of colors to use for the plots.")
    parser.add_argument("--positions_out_ORF", action="store", dest="positions_out_ORF", type=int, default=50\
                                             , help="The number of positions upstream of the start codon to include in the metagene vector. Default: 20.")
    parser.add_argument("--positions_in_ORF", action="store", dest="positions_in_ORF", type=int, default=200\
                                            , help="The number of positions downstream of the start codon to include in the metagene vector. Default: 100.")
    parser.add_argument("--output_formats", nargs="+", action="store", dest="output_formats", default=["interactive", "svg"]\
                                            , help="The output format of the plots (interactive, svg, pdf, jpg, png). Default: interactive, svg.")
    parser.add_argument("--include_plotly_js", action="store", dest="include_plotly_js", type=str, default="integrated",\
                                            help="The way to include the plotly.js library (integrated, local, online). Default: integrated.")
    args = parser.parse_args()

    alignment_file = args.alignment_file_path
    genome_length_dict = io.parse_genome_lengths(args.genome_file_path)
    mapping_methods = io.parse_mapping_methods(args.mapping_methods)
    read_length_list = io.parse_read_lengths(args.read_lengths)

    fig_list = []
    ir = IntervalReader(alignment_file)
    read_intervals_dict, total_counts_dict = ir.output()

    for normalization_method in args.normalization_methods:
        meta_dir = args.output_dir_path / alignment_file.stem / args.normalization_method
        meta_dir.mkdir(parents=True, exist_ok=True)

        for mapping_method in mapping_methods:
            start_codon_dict, stop_codon_dict, _ = ann.retrieve_annotation_positions(args.annotation_file_path, read_intervals_dict, total_counts_dict, genome_length_dict, args.filtering_methods,\
                                                                                        mapping_method, args.rpkm_threshold, args.neighboring_genes_distance, args.positions_out_ORF, args.positions_in_ORF)

            start_coverage_dict = mg.metagene_mapping_start(start_codon_dict, read_intervals_dict, args.positions_out_ORF, args.positions_in_ORF, mapping_method)
            stop_coverage_dict = mg.metagene_mapping_stop(stop_codon_dict, read_intervals_dict, args.positions_out_ORF, args.positions_in_ORF, mapping_method)

            if normalization_method == "cpm":
                start_coverage_dict = misc.normalize_coverage(start_coverage_dict, total_counts_dict)
                stop_coverage_dict = misc.normalize_coverage(stop_coverage_dict, total_counts_dict)

            tmp = create_metagene_figures(start_coverage_dict, read_length_list, meta_dir, mapping_method, args.normalization_method,\
                                                                args.positions_out_ORF, args.positions_in_ORF, args.color_list)
            fig_list.extend(tmp)

        io.write_plots_to_file(fig_list, args.output_formats, args.include_plotly_js, alignment_file.stem, meta_dir)

if __name__ == '__main__':
    main()
