#!/usr/bin/env python

import argparse
import html
from pathlib import Path

import numpy as np
import pandas as pd

import lib.io as io
import lib.misc as misc
import lib.annotation as ann
import lib.metagene as mg
import lib.plotting as plotting
import lib.psite as psite

from lib.alignment import IntervalReader


ALL_CONTIGS = "[all contigs]"


def summarize_candidate_counts(selection, mapping_method):
    """Count unique CDSs and first-failure exclusions, including empty input."""
    groups = [(ALL_CONTIGS, selection)]
    groups.extend(selection.groupby("contig", sort=True))
    rows = []
    for contig, candidates in groups:
        reasons = candidates["reason"]
        rows.append({
            "mapping_method": mapping_method,
            "contig": contig,
            "input_cds": len(candidates),
            "cohort_cds": int((reasons != "cohort").sum()),
            "retained_cds": int((candidates["status"] == "retained").sum()),
            **{
                f"excluded_{reason}": int((reasons == reason).sum())
                for reason in ("cohort", "overlap", "length", "rpkm", "boundary", "error")
            },
        })
    return pd.DataFrame(rows)


def summarize_candidate_support(
    selection, read_intervals, read_lengths, mapping_method,
    positions_out_ORF, positions_in_ORF, anchors=("start", "stop"),
):
    """Count actual contributions to each anchor, rather than interval overlaps.

    Endpoint profiles count only an endpoint inside the window. Global profiles
    count aligned nucleotides. A CDS is counted once in the pooled selected-length
    row even if several lengths contribute. Selection coordinates are 0-based
    inclusive, as in lib.annotation and lib.metagene.
    """
    lengths = sorted(set(read_lengths))
    selected_lengths = set(lengths)
    contigs = sorted(set(selection["contig"]))
    support = {}
    for contig in [ALL_CONTIGS, *contigs]:
        retained = selection[selection["status"] == "retained"]
        if contig != ALL_CONTIGS:
            retained = retained[retained["contig"] == contig]
        for anchor in anchors:
            for length in [*lengths, "all_selected"]:
                support[(contig, anchor, length)] = {
                    "mapping_method": mapping_method,
                    "contig": contig,
                    "anchor": anchor,
                    "read_length": length,
                    "retained_cds": len(retained),
                    "contributing_cds": 0,
                    "raw_count_contributions": 0,
                }

    retained = selection[selection["status"] == "retained"]
    for candidate in retained.itertuples(index=False):
        read_index = read_intervals.get((candidate.contig, candidate.strand))
        if read_index is None:
            continue
        windows = ann.metagene_window_bounds(
            candidate.start, candidate.end, candidate.strand,
            positions_out_ORF, positions_in_ORF,
        )
        for anchor in anchors:
            low, high = windows[anchor]
            contributions = {length: 0 for length in lengths}
            for interval in read_index.find((low, high)):
                length = interval[2]
                if length not in selected_lengths:
                    continue
                if mapping_method == "global":
                    amount = sum(
                        max(0, min(high, end) - max(low, start) + 1)
                        for start, end in misc.get_aligned_blocks(interval)
                    )
                else:
                    position = mg._point_position(interval, candidate.strand, mapping_method)
                    amount = int(position is not None and low <= position <= high)
                contributions[length] += amount

            for contig in (candidate.contig, ALL_CONTIGS):
                for length, amount in contributions.items():
                    support[(contig, anchor, length)]["contributing_cds"] += int(amount > 0)
                    support[(contig, anchor, length)]["raw_count_contributions"] += amount
                pooled = support[(contig, anchor, "all_selected")]
                pooled["contributing_cds"] += int(any(contributions.values()))
                pooled["raw_count_contributions"] += sum(contributions.values())
    return pd.DataFrame(support.values())


def candidate_report_html(counts, support, description):
    """Embed candidate numbers and their definitions alongside the profiles."""
    return (
        f"<p>{html.escape(description)}</p>"
        "<p>Counts refer to unique CDS intervals. Retained CDSs pass the filters; "
        "contributing CDSs supply at least one actual count within the displayed "
        "window at the selected read length. All-selected counts count each CDS "
        "once across the selected lengths. Exclusions record the first failed "
        "filter. Candidate coordinates in the TSV are zero-based and inclusive.</p>"
        '<p>Download <a href="../candidates.tsv">candidate identities and exclusions</a>, '
        '<a href="../candidate_counts.tsv">filter counts</a>, or '
        '<a href="../candidate_support.tsv">profile support counts</a>.</p>'
        '<h2>Candidate counts</h2><div class="table-wrap">'
        + counts.to_html(index=False, border=0)
        + '</div><h2>Profile support</h2><div class="table-wrap">'
        + support.to_html(index=False, border=0)
        + "</div>"
    )


def estimate_profile_offsets(df_start, coordinates, mapping_method):
    """Estimate P-site markers only for profiles anchored on a physical read end."""
    columns = df_start.columns[1:]
    if mapping_method not in psite.READ_ENDS:
        return {int(column): None for column in columns}

    read_end = psite.READ_ENDS[mapping_method]
    offsets = {}
    for column in columns:
        estimate = psite.estimate_offset(
            df_start[column].to_numpy(dtype=float),
            coordinates,
            read_end,
        )
        offsets[int(column)] = estimate.offset if estimate.is_significant else None
    return offsets


def prepare_metagene_dataframes(
    start_coverage_dict,
    stop_coverage_dict,
    read_length_list,
    positions_out_ORF,
    positions_in_ORF,
):
    """Select configured lengths and construct balanced start/stop frames."""
    start_coverage_dict = misc.retain_read_lengths(
        start_coverage_dict, read_length_list
    )
    stop_coverage_dict = misc.retain_read_lengths(
        stop_coverage_dict, read_length_list
    )
    start_coverage_dict, stop_coverage_dict = misc.equalize_dictionary_keys(
        start_coverage_dict,
        stop_coverage_dict,
        positions_out_ORF,
        positions_in_ORF,
        read_lengths=read_length_list,
    )

    df_start_dict = misc.create_data_frame(
        start_coverage_dict,
        positions_out_ORF,
        positions_in_ORF,
        "start",
        read_lengths=read_length_list,
    )
    df_stop_dict = misc.create_data_frame(
        stop_coverage_dict,
        positions_out_ORF,
        positions_in_ORF,
        "stop",
        read_lengths=read_length_list,
    )
    return df_start_dict, df_stop_dict


def estimate_coverage_offsets(
    start_coverage_dict,
    stop_coverage_dict,
    read_length_list,
    mapping_method,
    positions_out_ORF,
    positions_in_ORF,
):
    """Estimate markers from raw counts before presentation normalization."""
    df_start_dict, _ = prepare_metagene_dataframes(
        start_coverage_dict,
        stop_coverage_dict,
        read_length_list,
        positions_out_ORF,
        positions_in_ORF,
    )
    coordinates = np.arange(-positions_out_ORF, positions_in_ORF)
    return {
        chromosome: estimate_profile_offsets(
            df_start, coordinates, mapping_method
        )
        for chromosome, df_start in df_start_dict.items()
    }


def create_metagene_figures(
    start_coverage_dict,
    stop_coverage_dict,
    read_length_list,
    meta_dir,
    mapping_method,
    normalization_method,
    positions_out_ORF,
    positions_in_ORF,
    color_list,
    offsets_by_chromosome=None,
    candidate_counts=None,
    candidate_support=None,
    start_only=False,
):
    """Create figures for one mapping and presentation normalization."""
    df_start_dict, df_stop_dict = prepare_metagene_dataframes(
        start_coverage_dict,
        stop_coverage_dict,
        read_length_list,
        positions_out_ORF,
        positions_in_ORF,
    )

    window_size = positions_out_ORF + positions_in_ORF
    coordinates = np.arange(-positions_out_ORF, positions_in_ORF)
    if offsets_by_chromosome is None:
        offsets_by_chromosome = {
            chromosome: estimate_profile_offsets(
                df_start, coordinates, mapping_method
            )
            for chromosome, df_start in df_start_dict.items()
        }

    fig_list = []
    for chromosome in df_start_dict:
        df_start = df_start_dict[chromosome]
        df_stop = df_stop_dict[chromosome]

        if normalization_method == "window":
            df_start = misc.window_normalize_df(df_start, window_size)
            df_stop = misc.window_normalize_df(df_stop, window_size)

        # Only read-end profiles have a biologically defined P-site geometry.
        # Centered and global profiles still retain the keys expected by the
        # plotting functions, but deliberately carry no offset markers.
        offsets = offsets_by_chromosome[chromosome]

        subtitle = f"{mapping_method} mapping, {normalization_method} normalisation"
        panel_support = None
        if candidate_counts is not None and candidate_support is not None:
            support_contig = ALL_CONTIGS if chromosome == "no_evidence" else chromosome
            counts = candidate_counts[candidate_counts["contig"] == support_contig]
            support = candidate_support[
                (candidate_support["contig"] == support_contig)
                & (candidate_support["anchor"] == "start")
            ]
            retained_count = int(counts["retained_cds"].sum())
            pooled_count = int(support.loc[
                support["read_length"] == "all_selected", "contributing_cds"
            ].sum())
            subtitle += f"\n{retained_count} retained CDSs · {pooled_count} contributing CDSs at starts"
            panel_support = {
                int(row.read_length): {
                    "contributing_cds": row.contributing_cds,
                    "raw_count_contributions": row.raw_count_contributions,
                }
                for row in support.itertuples(index=False)
                if row.read_length != "all_selected"
            }
        if not start_only:
            fig = plotting.plot_metagene_heatmap(
                df_start,
                df_stop,
                read_length_list,
                chromosome,
                subtitle,
                offsets,
                read_end=mapping_method,
            )
            fig_list.append((chromosome, mapping_method, fig))

        value_label = {"raw": "Reads", "cpm": "CPM", "window": "Window-normalized counts"}[normalization_method]
        overlaid = plotting.plot_metagene_profiles(
            df_start,
            df_stop,
            read_length_list,
            f"{chromosome}: overlaid read lengths",
            subtitle,
            color_list=color_list,
            value_label=value_label,
            start_only=start_only,
        )
        if overlaid is not None:
            fig_list.append((f"{chromosome} (overlaid read lengths)", mapping_method, overlaid))

        profiles = plotting.plot_read_length_profiles(
            df_start,
            read_length_list,
            f"{chromosome}: start codon profiles",
            subtitle,
            offsets,
            read_end=mapping_method,
            color_list=color_list,
            max_panels=len(read_length_list),
            candidate_support=panel_support,
            value_label=value_label,
        )
        if profiles is not None:
            fig_list.append((f"{chromosome} (per read length)", mapping_method, profiles))

    io.create_excel_file(df_start_dict, meta_dir / f"{mapping_method}_readcounts_start.xlsx")
    if not start_only:
        io.create_excel_file(df_stop_dict, meta_dir / f"{mapping_method}_readcounts_stop.xlsx")

    return fig_list


def profile_cohort(args, read_intervals, total_counts, genome_lengths, max_length=None):
    """Profile the ordinary CDS set or an additional, length-decoupled sORF set."""
    is_sorf = max_length is not None
    root = args.output_dir_path / "sorfs" if is_sorf else args.output_dir_path
    root.mkdir(parents=True, exist_ok=True)
    methods = io.parse_mapping_methods(args.mapping_methods)
    lengths = io.parse_read_lengths(args.read_lengths)
    filters = [method for method in args.filtering_methods if method != "length"] if is_sorf else args.filtering_methods
    profile_data = {}
    selections, count_tables, support_tables = [], [], []
    for mapping_method in methods:
        starts, stops, selection = ann.retrieve_annotation_positions(
            args.annotation_file_path, read_intervals, total_counts, genome_lengths,
            filters, mapping_method, args.rpkm_threshold, args.neighboring_genes_distance,
            args.positions_out_ORF, args.positions_in_ORF, args.length_cutoff,
            cds_max_length=max_length, return_selection=True,
            required_anchors=("start",) if is_sorf else ("start", "stop"),
        )
        print(f"{'sORF' if is_sorf else 'General'} profile: {mapping_method}")
        start = mg.metagene_mapping_start(starts, read_intervals, args.positions_out_ORF, args.positions_in_ORF, mapping_method)
        stop = {} if is_sorf else mg.metagene_mapping_stop(stops, read_intervals, args.positions_out_ORF, args.positions_in_ORF, mapping_method)
        offsets = estimate_coverage_offsets(start, stop, lengths, mapping_method, args.positions_out_ORF, args.positions_in_ORF)
        counts = summarize_candidate_counts(selection, mapping_method)
        support = summarize_candidate_support(
            selection, read_intervals, lengths, mapping_method,
            args.positions_out_ORF, args.positions_in_ORF,
            anchors=("start",) if is_sorf else ("start", "stop"),
        )
        selections.append(selection.assign(mapping_method=mapping_method))
        count_tables.append(counts)
        support_tables.append(support)
        profile_data[mapping_method] = (start, stop, offsets, counts, support)

    candidates = pd.concat(selections, ignore_index=True)
    candidates = candidates[["mapping_method", *ann.SELECTION_COLUMNS]]
    counts = pd.concat(count_tables, ignore_index=True)
    support = pd.concat(support_tables, ignore_index=True)
    for filename, table in (("candidates.tsv", candidates), ("candidate_counts.tsv", counts), ("candidate_support.tsv", support)):
        table.to_csv(root / filename, sep="\t", index=False)

    description = (
        f"sORF group: reference-annotation CDSs with length < {max_length} nt. "
        "The ordinary minimum-length filter is disabled for this group. Fixed start-relative "
        "windows may extend beyond a short CDS into surrounding sequence."
        if is_sorf else "General group: reference-annotation CDSs passing the configured filters."
    )
    report = candidate_report_html(counts, support, description)
    for normalization_method in args.normalization_methods:
        normalization_method = normalization_method.lower()
        meta_dir = root / normalization_method
        meta_dir.mkdir(parents=True, exist_ok=True)
        figures = []
        for mapping_method, (start, stop, offsets, method_counts, method_support) in profile_data.items():
            # Presentation normalization must not alter later normalizations or
            # the raw counts used for support and P-site markers.
            presented_start = {contig: {length: values.copy() for length, values in coverage.items()} for contig, coverage in start.items()}
            presented_stop = {contig: {length: values.copy() for length, values in coverage.items()} for contig, coverage in stop.items()}
            if normalization_method == "cpm":
                presented_start = misc.normalize_coverage(presented_start, total_counts)
                presented_stop = misc.normalize_coverage(presented_stop, total_counts)
            figures.extend(create_metagene_figures(
                presented_start, presented_stop, lengths, meta_dir, mapping_method,
                normalization_method, args.positions_out_ORF, args.positions_in_ORF,
                args.color_list, offsets_by_chromosome=offsets,
                candidate_counts=method_counts, candidate_support=method_support,
                start_only=is_sorf,
            ))
        io.write_plots_to_file(
            figures, args.output_formats, args.include_plotly_js,
            args.alignment_file_path.stem + (" (sORFs)" if is_sorf else ""), meta_dir,
            report_html=report,
            report_subtitle=description + f" Normalization: {normalization_method}. "
                + ("Start-codon read-end profiles." if is_sorf else "Heatmap colours show enrichment relative to each length's background; line profiles use the stated normalization."),
        )


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
    parser.add_argument("--sorf_max_length", type=int, default=0,
                        help="Add start-codon profiles for CDSs shorter than this many nt, without the ordinary minimum-length filter; 0 disables.")
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
    parser.add_argument("--output_formats", nargs="+", action="store", dest="output_formats", default=["interactive", "svg"],
                                            choices=["interactive", "svg", "pdf", "jpg", "png"]\
                                            , help="The output format of the plots (interactive, svg, pdf, jpg, png). Default: interactive, svg.")
    parser.add_argument("--include_plotly_js", action="store", dest="include_plotly_js", type=str, default="integrated",\
                                            help="The way to include the plotly.js library (integrated, local, online). Default: integrated.")
    args = parser.parse_args()
    if args.sorf_max_length < 0:
        parser.error("--sorf_max_length must be nonnegative")

    alignment_file = args.alignment_file_path
    genome_length_dict = io.parse_genome_lengths(args.genome_file_path)
    ir = IntervalReader(alignment_file)
    read_intervals_dict, total_counts_dict = ir.output()
    profile_cohort(args, read_intervals_dict, total_counts_dict, genome_length_dict)
    if args.sorf_max_length:
        profile_cohort(args, read_intervals_dict, total_counts_dict, genome_length_dict, args.sorf_max_length)

if __name__ == '__main__':
    main()
