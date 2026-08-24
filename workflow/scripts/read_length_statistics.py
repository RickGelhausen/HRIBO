#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

import lib.theme as theme


from lib.alignment import LengthCounter
import lib.io as io

def parse_alignment_files(alignment_dir_path):
    """
    Read alignment files from directory.
    """

    files = [entry for entry in Path(alignment_dir_path).glob("*.bam") if entry.is_file()]
    files.extend([entry for entry in Path(alignment_dir_path).glob("*.sam") if entry.is_file()])

    return sorted(files, key=lambda x: str(x).lower())


def data_equalizer(read_length_counts, possible_alignment_files, possible_read_lengths):
    """
    Ensure that for each existing chromosome all alignment files and all read length exist,
    otherwise add 0 counts to make processing easier.
    """

    for chrom in read_length_counts:
        for possible_aligment_file in possible_alignment_files:
            if possible_aligment_file.stem not in read_length_counts[chrom]:
                read_length_counts[chrom][possible_aligment_file.stem] = {}
                for possible_read_length in possible_read_lengths[chrom]:
                    read_length_counts[chrom][possible_aligment_file.stem][possible_read_length] = 0
            else:
                for possible_read_length in possible_read_lengths[chrom]:
                    if possible_read_length not in read_length_counts[chrom][possible_aligment_file.stem]:
                        read_length_counts[chrom][possible_aligment_file.stem][possible_read_length] = 0

    return read_length_counts

def create_dataframes(read_length_counts):
    """
    Create dataframes from read length counts.
    """

    dataframes = {}
    for chrom in read_length_counts:
        for alignment_file in read_length_counts[chrom]:
            if chrom not in dataframes:
                dataframes[chrom] = pd.DataFrame()
                dataframes[chrom]["read_lengths"] = [str(x) for x in sorted(read_length_counts[chrom][alignment_file])]
                dataframes[chrom][alignment_file] = [int(read_length_counts[chrom][alignment_file][x]) for x in sorted(read_length_counts[chrom][alignment_file])]
            else:
                dataframes[chrom][alignment_file] = [int(read_length_counts[chrom][alignment_file][x]) for x in sorted(read_length_counts[chrom][alignment_file])]

    return dataframes


def create_dataframes_fractions(read_length_counts):
    """
    Create dataframes from read length counts.
    """

    dataframes = {}
    for chrom in read_length_counts:
        for alignment_file in read_length_counts[chrom]:

            read_counts = [int(read_length_counts[chrom][alignment_file][x]) for x in sorted(read_length_counts[chrom][alignment_file])]
            read_counts = [x / sum(read_counts) if sum(read_counts) > 0 else 0 for x in read_counts]
            if chrom not in dataframes:
                dataframes[chrom] = pd.DataFrame()
                dataframes[chrom]["read_lengths"] = [str(x) for x in sorted(read_length_counts[chrom][alignment_file])]
                dataframes[chrom][alignment_file] = read_counts
            else:
                dataframes[chrom][alignment_file] = read_counts

    return dataframes

def create_excel_output(dataframes, output_file):
    """
    Create excel output.
    """

    out_df = {}
    for chrom in dataframes:
        out_df[f"{chrom}"] = dataframes[chrom]

    io.excel_writer(output_file, out_df)

def plot_length_fractions(dataframes, color_list):
    """One read length distribution per sequence, all libraries overlaid.

    Line rather than area: the question is which read lengths a library is
    enriched for, and overlaying a handful of libraries answers it directly.
    """
    figure_dict = {}

    for chrom in dataframes:
        fig = go.Figure()
        cur_df = dataframes[chrom]
        labels = cur_df["read_lengths"]

        for index, file in enumerate(cur_df.columns[1:]):
            data = cur_df[file]
            data = data / data.sum()
            fig.add_trace(
                go.Scatter(
                    x=labels,
                    y=data,
                    name=file,
                    mode="lines",
                    line=dict(color=theme.series_color(index), width=2),
                    hovertemplate=f"{file}<br>read length %{{x}} nt<br>%{{y:.1%}} of reads<extra></extra>",
                )
            )

        fig.update_xaxes(title_text="Read length (nt)", type="category")
        fig.update_yaxes(title_text="Fraction of mapped reads", tickformat=".0%")

        theme.apply(
            fig,
            chrom,
            "Share of mapped reads at each length, per library",
        )
        fig.update_layout(
            height=460,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        )
        figure_dict[chrom] = fig

    return figure_dict


def create_html_output(figure_dict, output_path):
    """Render every figure into one page, with plotly.js included once."""
    parts = []
    for index, chrom in enumerate(figure_dict):
        parts.append('<div class="plot">')
        parts.append(
            figure_dict[chrom].to_html(
                full_html=False,
                include_plotlyjs=(index == 0),
                default_width="100%",
                config={"displaylogo": False, "responsive": True},
            )
        )
        parts.append("</div>")

    with open(output_path / "read_length_fractions.html", "w") as f:
        f.write(
            theme.page(
                "Read length distribution",
                "Fraction of mapped reads at each read length, for every library.",
                "\n".join(parts),
            )
        )


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Create plot with fractional readcounts", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-a", "--alignment_files", nargs="+", dest="alignment_files", help="Path to bam files.")
    parser.add_argument("-r", "--read_lengths", action="store", dest="read_lengths", type=str, default="10-40"\
                                              , help="The read lengths to be considered for the metagene-profiling.")
    parser.add_argument("-o", "--output_dir_path", action="store", dest="output_dir_path", type=Path, required=True\
                                                 , help="The path to the output directory.")
    parser.add_argument("--color_list", nargs="+", dest="color_list", required=False,\
                                        default=None, help="List of colors to use for the plots. Default contains 7 colorblind-friendly colors.")
    args = parser.parse_args()

    alignment_files = [Path(x) for x in args.alignment_files]
    output_dir_path = args.output_dir_path
    output_dir_path.mkdir(parents=True, exist_ok=True)
    read_length_list = io.parse_read_lengths(args.read_lengths)

    read_length_counts = {}
    possible_read_length_dict = {}
    for alignment_file in alignment_files:
        lc = LengthCounter(alignment_file, read_length_list)
        length_count_dict = lc.output()
        for chrom in length_count_dict.keys():
            if chrom not in read_length_counts:
                read_length_counts[chrom] = {}
            if chrom not in possible_read_length_dict:
                possible_read_length_dict[chrom] = set()

            for read_length in length_count_dict[chrom]:
                if read_length not in possible_read_length_dict[chrom]:
                    possible_read_length_dict[chrom].add(read_length)

            read_length_counts[chrom][alignment_file.stem] = length_count_dict[chrom]


    read_length_counts = data_equalizer(read_length_counts, alignment_files, possible_read_length_dict)

    dataframes = create_dataframes(read_length_counts)
    dataframes_fraction = create_dataframes_fractions(read_length_counts)
    create_excel_output(dataframes, output_dir_path / "read_length_counts.xlsx")
    create_excel_output(dataframes_fraction, output_dir_path / "read_length_fractions.xlsx")
    figure_dict = plot_length_fractions(dataframes, args.color_list)
    create_html_output(figure_dict, output_dir_path)

if __name__ == '__main__':
    main()