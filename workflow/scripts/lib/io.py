"""
Contains scripts used for parsing input data and writing output data.
Author: Rick Gelhausen
"""

import sys
from pathlib import Path

import pandas as pd

from lib import theme

def parse_alignment_files(alignment_dir_path):
    """
    Read alignment files from directory.
    """

    files = [entry for entry in Path(alignment_dir_path).glob("*.bam") if entry.is_file()]
    files.extend([entry for entry in Path(alignment_dir_path).glob("*.sam") if entry.is_file()])

    return sorted(files, key=lambda x: str(x).lower())

def parse_read_lengths(read_lengths):
    """
    Parse the read length input into a continuous list form.
    """

    parts = read_lengths.split(",")
    read_lengths = set()
    for part in parts:
        if "-" in part:
            interval = part.split("-")
            if int(interval[0]) < int(interval[1]):
                i1, i2 = int(interval[0]), int(interval[1])
            else:
                i1, i2 = int(interval[1]), int(interval[0])

            for i in range(i1, i2+1):
                read_lengths.add(i)
        else:
            read_lengths.add(int(part))

    return sorted(read_lengths)

def parse_genome_lengths(genome_file_path):
    """
    Read  lengths from genome file.
    """

    genome_length_dict = {}
    with open(genome_file_path, "r") as genome_file:
        for line in genome_file:
            if line[0] == ">":
                chromosome = line[1:].split(" ")[0].strip()
                genome_length_dict[chromosome] = 0
            else:
                genome_length_dict[chromosome] += len(line.strip())

    return genome_length_dict

def parse_mapping_methods(mapping_methods):
    """
    Parse the mapping methods into a list.
    """

    for method in mapping_methods:
        if method not in ["fiveprime", "threeprime", "global", "centered"]:
            sys.exit(f"Error: mapping method {method} not recognized. Please use one of the following: fiveprime, threeprime, global, centred.")

    return mapping_methods

def create_excel_file(in_df, output_file):
    """
    Create a sheet for the excel file containing the metagene profiling read counts.
    """

    for chromosome in in_df:
        df = in_df[chromosome]
        in_df[chromosome]["sum"] = df[df.columns[1:]].sum(numeric_only=True, axis=1)

    excel_writer(output_file, in_df)

def excel_writer(output_path, data_frames):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  []
    writer = pd.ExcelWriter(output_path, engine='xlsxwriter')
    for sheetname, df in data_frames.items():
        df.to_excel(writer, sheet_name=sheetname, index=False)
        worksheet = writer.sheets[sheetname]
        worksheet.freeze_panes(1, 0)
        for idx, col in enumerate(df):
            series = df[col]
            if col in header_only:
                max_len = len(str(series.name)) + 2
            else:
                max_len = max(( series.astype(str).str.len().max(), len(str(series.name)) )) + 1
            worksheet.set_column(idx, idx, max_len)
    writer.close()

def write_plots_to_file(fig_list, output_format, include_plotly_js, alignment_file_name, meta_dir, fig_width=1400, fig_height=600):
    """
    Write plots to requested file formats
    """

    if "png" in output_format:
        for chromosome, mapping_method, fig in fig_list:
            fig.write_image(f"{meta_dir}/{chromosome}_{mapping_method}.png", width=fig_width, height=fig_height)

    if "jpg" in output_format:
        for chromosome, mapping_method, fig in fig_list:
            fig.write_image(f"{meta_dir}/{chromosome}_{mapping_method}.jpg", width=fig_width, height=fig_height)

    if "svg" in output_format:
        for chromosome, mapping_method, fig in fig_list:
            fig.write_image(f"{meta_dir}/{chromosome}_{mapping_method}.svg", width=fig_width, height=fig_height)

    if "interactive" in output_format:
        create_interactive_html(fig_list, alignment_file_name, f"{meta_dir}/interactive_metagene_profiling.html", include_plotly_js)

def create_interactive_html(fig_list, alignment_file_name, output_file, include_plotly_js):
    """Render every figure into one standalone page.

    plotly.js is emitted once for the whole page rather than once per figure,
    which is what made the previous reports grow by several megabytes for each
    additional plot.
    """

    mapping_labels = {
        "global": "Global mapping",
        "threeprime": "3' mapping",
        "fiveprime": "5' mapping",
        "centered": "Centered mapping",
    }

    js_mode = {"integrated": True, "online": "cdn", "local": "directory"}.get(include_plotly_js, True)

    parts = []
    seen_headings = []
    for index, (name, mapping_method, fig) in enumerate(fig_list):
        heading = mapping_labels.get(mapping_method, mapping_method)
        if heading not in seen_headings:
            parts.append(f"<h2>{heading}</h2>")
            seen_headings.append(heading)

        parts.append(f"<h3>{name}</h3>")
        parts.append('<div class="plot">')
        parts.append(
            fig.to_html(
                full_html=False,
                include_plotlyjs=(js_mode if index == 0 else False),
                default_width="100%",
                config={"displaylogo": False, "responsive": True},
            )
        )
        parts.append("</div>")

    with open(output_file, "w") as handle:
        handle.write(
            theme.page(
                f"Metagene profiling: {alignment_file_name}",
                "Enrichment is shown relative to each read length's own background, "
                "so that a sparse read length stays legible beside a deep one.",
                "\n".join(parts),
            )
        )
