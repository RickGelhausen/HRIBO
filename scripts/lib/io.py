"""
Contains scripts used for parsing input data and writing output data.
Author: Rick Gelhausen
"""

import sys
import pandas as pd
from pathlib import Path

INTRO_HTML = \
"""
<!DOCTYPE html>
<html>
    <head>
        <style>
            h1 { text-align: center; font-size: 3.5em; }
            h2 { text-align: center; font-size: 2em; }
            h3 { text-align: center; font-size: 1.5em; }
            h4 { text-align: left; font-size: 1.5em; }
            .headerline { height: 2px; background-color: black; }
            .seperator { height: 1px; width: 90%;}
            .page { margin: 50px; text-align: center; }
            .plotly-graph-div { margin: 0 auto; }
            .toggleBtn { margin: 20px; }
        </style>
        <title>Metagene plots</title>
    </head>
    <body>
        <div class="page">
"""

OUTRO_HTML = \
"""
        </div>
    </body>
</html>
"""

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
            read_lengths.add(part)

    return [int(i) for i in sorted(list(read_lengths))]

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
    writer.save()

def write_plots_to_file(fig_list, output_format, include_plotly_js, alignment_file_name, meta_dir):
    """
    Write plots to requested file formats
    """

    if "png" in output_format:
        for chromosome, mapping_method, fig in fig_list:
            fig.write_image(f"{meta_dir}/{chromosome}_{mapping_method}.png")

    if "jpg" in output_format:
        for chromosome, mapping_method, fig in fig_list:
            fig.write_image(f"{meta_dir}/{chromosome}_{mapping_method}.jpg")

    if "svg" in output_format:
        for chromosome, mapping_method, fig in fig_list:
            fig.write_image(f"{meta_dir}/{chromosome}_{mapping_method}.svg")

    if "interactive" in output_format:
        create_interactive_html(fig_list, alignment_file_name, f"{meta_dir}/interactive_metagene_profiling.html", include_plotly_js)

def create_interactive_html(fig_list, alignment_file_name, output_file, include_plotly_js):
    """
    Takes a list of figures and creates an output HTML form.
    """

    in_header = []
    html_string = INTRO_HTML + "\n"
    html_string += f"<h1>{alignment_file_name}</h1>\n"
    for idx, (chromosome, mapping_method, fig) in enumerate(fig_list):
        if mapping_method == "global":
            mapping_out = "Global Mapping"
        elif mapping_method == "threeprime":
            mapping_out = "3' Mapping"
        elif mapping_method == "fiveprime":
            mapping_out = "5' Mapping"
        elif mapping_method == "centered":
            mapping_out = "Centered Mapping"
        else:
            mapping_out = "Unknown"

        if mapping_out not in in_header:
            html_string += "<hr class=headerline>\n"
            html_string += f"<h2>{mapping_out}</h2>\n"
            in_header.append(mapping_out)

        if include_plotly_js == "integrated":
            if idx == 0:
                html_string += "<div class=plot>\n"
                html_string += fig.to_html(full_html=False, default_height="600px", default_width="80%")
                html_string += "</div>\n"
            else:
                html_string += "<div class=plot>\n"
                html_string += fig.to_html(full_html=False, include_plotlyjs=False, default_height="600px", default_width="80%")
                html_string += "</div>\n"
        elif include_plotly_js == "online":
            if idx == 0:
                html_string += "<div class=plot>\n"
                html_string += fig.to_html(full_html=False, include_plotlyjs="cdn", default_height="600px", default_width="80%")
                html_string += "</div>\n"
            else:
                html_string += "<div class=plot>\n"
                html_string += fig.to_html(full_html=False, include_plotlyjs=False, default_height="600px", default_width="80%")
                html_string += "</div>\n"
        elif include_plotly_js == "local":
            if idx == 0:
                html_string += "<div class=plot>\n"
                html_string += fig.to_html(full_html=False, include_plotlyjs="directory", default_height="600px", default_width="80%")
                html_string += "</div>\n"
            else:
                html_string += "<div class=plot>\n"
                html_string += fig.to_html(full_html=False, include_plotlyjs=False, default_height="600px", default_width="80%")
                html_string += "</div>\n"



    html_string += OUTRO_HTML
    with open(output_file, "w") as f:
        f.write(html_string)

