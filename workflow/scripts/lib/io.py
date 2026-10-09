"""
Contains scripts used for parsing input data and writing output data.
Author: Rick Gelhausen
"""

import hashlib
import html
import json
import re
import sys
from pathlib import Path
from urllib.parse import quote

import pandas as pd
from plotly.offline import get_plotlyjs

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

def write_plots_to_file(fig_list, output_format, include_plotly_js, alignment_file_name, meta_dir, fig_width=1400, fig_height=None, report_html="", report_subtitle=None, contig_reports=None):
    """
    Write plots to requested file formats
    """

    for image_format in ("png", "jpg", "svg", "pdf"):
        if image_format in output_format:
            for chromosome, mapping_method, fig in fig_list:
                output_path = Path(meta_dir) / f"{chromosome}_{mapping_method}.{image_format}"
                height = fig_height if fig_height is not None else (fig.layout.height or 600)
                fig.write_image(
                    str(output_path), width=fig_width, height=height
                )

    if "interactive" in output_format:
        create_interactive_html(fig_list, alignment_file_name, f"{meta_dir}/interactive_metagene_profiling.html", include_plotly_js, report_html, report_subtitle, contig_reports)

def metagene_contig(name, fig):
    """Read the explicit contig identity, with a fallback for older callers."""
    metadata = getattr(getattr(fig, "layout", None), "meta", None)
    if isinstance(metadata, dict) and "metagene_contig" in metadata:
        return str(metadata["metagene_contig"])
    for suffix in (" (overlaid read lengths)", " (per read length)"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def _metagene_figures_html(fig_list, js_mode):
    """Render one contig's figures, including Plotly once for the page."""
    mapping_labels = {
        "global": "Global mapping",
        "threeprime": "3' mapping",
        "fiveprime": "5' mapping",
        "centered": "Centered mapping",
    }

    parts = []
    seen_headings = []
    for index, (name, mapping_method, fig) in enumerate(fig_list):
        heading = mapping_labels.get(mapping_method, mapping_method)
        if heading not in seen_headings:
            parts.append(f"<h2>{heading}</h2>")
            seen_headings.append(heading)

        parts.append(f"<h3>{html.escape(name)}</h3>")
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

    return "\n".join(parts)


def create_interactive_html(fig_list, alignment_file_name, output_file, include_plotly_js, report_html="", report_subtitle=None, contig_reports=None):
    """Use an index and separate contig pages for multi-contig profiles.

    A single-contig report retains its existing inline layout. Contig pages
    stay beside the index so candidate-table download links remain valid.
    """
    output_file = Path(output_file)
    groups = {}
    for item in fig_list:
        name, _, fig = item
        groups.setdefault(metagene_contig(name, fig), []).append(item)

    js_mode = {"integrated": True, "online": "cdn", "local": "directory"}.get(include_plotly_js, True)
    if fig_list and js_mode == "directory":
        # to_html inserts the reference but, unlike write_html, does not copy JS.
        (output_file.parent / "plotly.min.js").write_text(get_plotlyjs())

    title = f"Metagene profiling: {alignment_file_name}"
    subtitle = report_subtitle if report_subtitle is not None else (
        "Enrichment is shown relative to each read length's own background, "
        "so that a sparse read length stays legible beside a deep one."
    )
    pages = {}
    if len(groups) > 1:
        for contig in sorted(groups):
            slug = re.sub(r"[^A-Za-z0-9_.-]", "_", contig)[:80] or "contig"
            digest = hashlib.sha256(contig.encode()).hexdigest()[:10]
            pages[contig] = f"{output_file.stem}_contig_{slug}-{digest}.html"

        def navigation(current=None):
            links = []
            for contig, filename in pages.items():
                label = "No profile evidence" if contig == "no_evidence" else contig
                attribute = ' aria-current="page"' if contig == current else ""
                links.append(
                    f'<li><a href="{html.escape(quote(filename))}"{attribute}>'
                    f'{html.escape(label)}</a></li>'
                )
            return '<nav aria-label="Contig profiles"><h2>Contig profiles</h2><ul>' + "".join(links) + "</ul></nav>"

        for contig, filename in pages.items():
            specific_report = (contig_reports or {}).get(contig, "")
            body = (
                f'<p><a href="{html.escape(quote(output_file.name))}">← Report index</a></p>'
                + navigation(contig)
                + _metagene_figures_html(groups[contig], js_mode)
                + specific_report
            )
            label = "No profile evidence" if contig == "no_evidence" else contig
            (output_file.parent / filename).write_text(
                theme.page(html.escape(f"{title} · {label}"), html.escape(subtitle), body)
            )
        body = (
            "<p>Choose a contig to open its profiles. Each report contains its "
            "configured mapping methods and plot views.</p>"
            + navigation()
            + report_html
        )
    else:
        body = report_html + _metagene_figures_html(fig_list, js_mode)
    output_file.write_text(theme.page(html.escape(title), html.escape(subtitle), body))

    # Only remove files recorded as this writer's generated pages on a prior run.
    manifest = output_file.with_suffix(".contigs.json")
    previous = json.loads(manifest.read_text()).get("pages", {}) if manifest.exists() else {}
    for filename in set(previous.values()) - set(pages.values()):
        if (Path(filename).name == filename
                and filename.startswith(f"{output_file.stem}_contig_")
                and filename.endswith(".html")):
            (output_file.parent / filename).unlink(missing_ok=True)
    manifest.write_text(json.dumps({"pages": pages}, indent=2) + "\n")
