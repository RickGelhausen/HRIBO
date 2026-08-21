#!/usr/bin/env python
"""
Recommend a translation initiation site caller setup for one library.

Given an alignment and an annotation, this works out which read lengths carry a
usable initiation signal, how far each of them sits from the P-site, and whether
the evidence is strong enough to act on. The result is written as a readable
HTML report, a machine-readable JSON document, and a TSV of the per-read-length
evidence.

The recommendation is deliberately allowed to be "none". Bacterial Ribo-seq
often lacks the periodicity that makes this easy in eukaryotes, and a confident
offset derived from noise is worse than no offset at all.

Author: Rick Gelhausen
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

import lib.annotation as ann
import lib.io as io
import lib.metagene as mg
import lib.misc as misc
import lib.plotting as plotting
import lib.psite as psite
import lib.theme as theme
from lib.alignment import IntervalReader


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Recommend read lengths and P-site offsets for a TIS caller.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-b", "--alignment_file_path", type=Path, required=True,
                        help="Alignment (.bam) file, named <METHOD>-<CONDITION>-<REPLICATE>.bam.")
    parser.add_argument("-a", "--annotation_file_path", type=Path, required=True,
                        help="Annotation file in GFF3 or GTF format.")
    parser.add_argument("-g", "--genome_file_path", type=Path, required=True,
                        help="Genome FASTA, used for the sequence lengths.")
    parser.add_argument("-o", "--output_dir_path", type=Path, required=True,
                        help="Directory for the report, JSON and TSV.")
    parser.add_argument("-r", "--read_lengths", type=str, default="22-40",
                        help="Read lengths to consider.")
    parser.add_argument("--mapping_method", type=str, default="fiveprime",
                        choices=["fiveprime", "threeprime"],
                        help="Read end the offsets are measured from.")
    parser.add_argument("--positions_out_ORF", type=int, default=100,
                        help="Nucleotides upstream of the start codon to profile.")
    parser.add_argument("--positions_in_ORF", type=int, default=150,
                        help="Nucleotides inside the ORF to profile.")
    parser.add_argument("--filtering_methods", nargs="+", default=["overlap", "rpkm", "length"],
                        help="Annotation filters applied before profiling.")
    parser.add_argument("--neighboring_genes_distance", type=int, default=50,
                        help="Distance used when removing overlapping genes.")
    parser.add_argument("--rpkm_threshold", type=float, default=10.0,
                        help="Minimum RPKM for a gene to contribute.")
    parser.add_argument("--include_plotly_js", type=str, default="integrated",
                        choices=["integrated", "online", "local"],
                        help="How the report references plotly.js. 'integrated' is self "
                             "contained but adds several megabytes per report.")
    return parser.parse_args()


def build_profiles(args):
    """Metagene start and stop profiles per read length, plus per-length totals."""
    genome_lengths = io.parse_genome_lengths(args.genome_file_path)
    read_lengths = io.parse_read_lengths(args.read_lengths)

    reader = IntervalReader(args.alignment_file_path)
    read_intervals, total_counts = reader.output()

    start_codons, stop_codons = ann.retrieve_annotation_positions(
        args.annotation_file_path,
        read_intervals,
        total_counts,
        genome_lengths,
        args.filtering_methods,
        args.mapping_method,
        args.rpkm_threshold,
        args.neighboring_genes_distance,
        args.positions_out_ORF,
        args.positions_in_ORF,
    )

    start_coverage = mg.metagene_mapping_start(
        start_codons, read_intervals, args.positions_out_ORF, args.positions_in_ORF, args.mapping_method
    )
    stop_coverage = mg.metagene_mapping_stop(
        stop_codons, read_intervals, args.positions_out_ORF, args.positions_in_ORF, args.mapping_method
    )
    start_coverage, stop_coverage = misc.equalize_dictionary_keys(
        start_coverage, stop_coverage, args.positions_out_ORF, args.positions_in_ORF
    )

    # Sum across sequences: the offset is a property of the protocol, not of a
    # particular replicon, and pooling gives the estimate more to work with.
    start_profiles = _sum_over_chromosomes(start_coverage, read_lengths)
    stop_profiles = _sum_over_chromosomes(stop_coverage, read_lengths)

    totals = {length: int(profile.sum()) for length, profile in start_profiles.items()}
    return start_profiles, stop_profiles, totals


def _sum_over_chromosomes(coverage, wanted_lengths):
    pooled: dict[int, np.ndarray] = {}
    for chromosome in coverage:
        for read_length, values in coverage[chromosome].items():
            if int(read_length) not in wanted_lengths:
                continue
            values = np.asarray(values, dtype=float)
            key = int(read_length)
            pooled[key] = values.copy() if key not in pooled else pooled[key] + values
    return pooled


def to_dataframe(profiles, coordinates):
    import pandas as pd

    frame = pd.DataFrame({"coordinates": coordinates})
    for read_length in sorted(profiles):
        frame[str(read_length)] = profiles[read_length]
    return frame


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------


def orfbounder_config(recommendation, mapping_method):
    """A config block that can be pasted straight into an ORFBounder run."""
    if not recommendation.has_recommendation:
        return "# No usable initiation signal was found; no setup is suggested."
    offsets = "\n".join(
        f"    {length}: {offset}" for length, offset in recommendation.offset_table()
    )
    return (
        "# Suggested ORFBounder settings, derived from this library.\n"
        f"readLengths: [{', '.join(str(l) for l in recommendation.read_lengths)}]\n"
        f"mappingMethod: \"{mapping_method}\"\n"
        "psiteOffsets:\n"
        f"{offsets}\n"
    )


def scores_to_tsv(scores, path):
    header = [
        "read_length", "total_reads", "abundance", "offset", "peak_height",
        "background", "sharpness", "z_score", "frame_0", "frame_1", "frame_2",
        "frame_bias", "periodicity", "usable", "reasons",
    ]
    with open(path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for score in scores:
            handle.write("\t".join([
                str(score.read_length),
                str(score.total_reads),
                f"{score.abundance:.6f}",
                "" if score.offset is None else str(score.offset),
                f"{score.peak_height:.2f}",
                f"{score.background:.2f}",
                f"{score.sharpness:.3f}",
                f"{score.z_score:.3f}",
                *[f"{value:.4f}" for value in score.frame_fractions],
                f"{score.frame_bias:.4f}",
                f"{score.periodicity:.4f}",
                "yes" if score.usable else "no",
                "; ".join(score.reasons),
            ]) + "\n")


def render_report(library, recommendation, scores, figures, mapping_method, path, include_plotly_js="integrated"):
    """The human-facing report: verdict first, then the evidence behind it."""
    colour = theme.CONFIDENCE_COLORS.get(recommendation.confidence, theme.INK_MUTED)

    parts = ['<div class="verdict">']
    if recommendation.has_recommendation:
        lengths = ", ".join(str(length) for length in recommendation.read_lengths)
        parts.append(
            f"<strong>Use read lengths {lengths}</strong>"
            f'<span class="pill" style="background:{colour}">{recommendation.confidence} confidence</span>'
        )
        parts.append(
            f"<p>The pooled initiation peak is {recommendation.sharpness:.1f} times the "
            f"upstream background, built from {recommendation.covered_fraction:.0%} of the "
            f"library. The dominant reading frame holds {recommendation.frame_bias:.0%} of "
            "P-sites.</p>"
        )
    else:
        parts.append(
            "<strong>No TIS caller setup is recommended for this library</strong>"
            f'<span class="pill" style="background:{colour}">no recommendation</span>'
        )
    parts.append("</div>")

    for warning in recommendation.warnings:
        parts.append(f'<p class="warn">{warning}</p>')

    parts.append("<h2>Suggested configuration</h2>")
    parts.append(f"<pre><code>{orfbounder_config(recommendation, mapping_method)}</code></pre>")

    if recommendation.rationale:
        parts.append("<h2>How this was chosen</h2><ul>")
        parts.extend(f"<li>{line}</li>" for line in recommendation.rationale)
        parts.append("</ul>")

    parts.append("<h2>Evidence per read length</h2>")
    parts.append(_scores_table(scores))

    parts.append("<h2>Figures</h2>")
    js_mode = {"integrated": True, "online": "cdn", "local": "directory"}.get(include_plotly_js, True)
    for index, (heading, figure) in enumerate(figures):
        if figure is None:
            continue
        parts.append(f"<h3>{heading}</h3>")
        parts.append('<div class="plot">')
        parts.append(
            figure.to_html(
                full_html=False,
                include_plotlyjs=(js_mode if index == 0 else False),
                default_width="100%",
                config={"displaylogo": False, "responsive": True},
            )
        )
        parts.append("</div>")

    Path(path).write_text(
        theme.page(
            f"TIS caller advice: {library}",
            f"{mapping_method} mapping. Offsets are the distance from the mapped read end "
            "to the P-site.",
            "\n".join(parts),
        )
    )


def _scores_table(scores):
    rows = [
        "<div class='table-wrap'><table><thead><tr>"
        "<th>Read length</th><th>Share of reads</th><th>Offset</th>"
        "<th>Peak vs background</th><th>Dominant frame</th><th>Periodicity</th>"
        "<th>Usable</th><th>Notes</th></tr></thead><tbody>"
    ]
    for score in scores:
        frame = "-" if not any(score.frame_fractions) else (
            f"{int(np.argmax(score.frame_fractions))} ({max(score.frame_fractions):.0%})"
        )
        rows.append(
            "<tr>"
            f"<td class='num'>{score.read_length}</td>"
            f"<td class='num'>{score.abundance:.1%}</td>"
            f"<td class='num'>{'-' if score.offset is None else score.offset}</td>"
            f"<td class='num'>{score.sharpness:.1f}x</td>"
            f"<td class='num'>{frame}</td>"
            f"<td class='num'>{score.periodicity:.0%}</td>"
            f"<td>{'yes' if score.usable else 'no'}</td>"
            f"<td>{'; '.join(score.reasons)}</td>"
            "</tr>"
        )
    rows.append("</tbody></table></div>")
    return "\n".join(rows)


def main():
    args = parse_arguments()
    args.output_dir_path.mkdir(parents=True, exist_ok=True)
    library = args.alignment_file_path.stem

    coordinates = np.arange(-args.positions_out_ORF, args.positions_in_ORF)
    start_profiles, stop_profiles, totals = build_profiles(args)

    if not start_profiles:
        scores, recommendation = [], psite._no_recommendation([])
    else:
        scores = psite.score_read_lengths(start_profiles, coordinates, totals)
        recommendation = psite.recommend_read_lengths(scores, start_profiles, coordinates)

    read_lengths = sorted(start_profiles)
    df_start = to_dataframe(start_profiles, coordinates)
    df_stop = to_dataframe(stop_profiles, coordinates)
    significant_offsets = {s.read_length: (s.offset if s.usable else None) for s in scores}

    figures = [
        ("Read length against position", plotting.plot_metagene_heatmap(
            df_start, df_stop, read_lengths, "Metagene profile",
            f"{args.mapping_method} mapping, enrichment over each read length's own background",
            significant_offsets)),
        ("Library composition", plotting.plot_read_length_distribution(
            scores, "Read length distribution",
            "Which read lengths carry a usable initiation signal")),
        ("Reading frame", plotting.plot_frame_composition(
            scores, "Reading frame composition",
            "Share of P-sites per frame, after offset correction")),
    ]
    if recommendation.has_recommendation:
        pooled = psite.pool_profiles(start_profiles, recommendation.offsets, recommendation.read_lengths)
        figures.append(("Pooled, offset corrected", plotting.plot_pooled_profile(
            pooled, coordinates, "Pooled P-site profile",
            f"read lengths {', '.join(str(l) for l in recommendation.read_lengths)}")))

    render_report(library, recommendation, scores, figures, args.mapping_method,
                  args.output_dir_path / "tis_recommendation.html", args.include_plotly_js)
    scores_to_tsv(scores, args.output_dir_path / "read_length_evidence.tsv")

    payload = {
        "library": library,
        "mapping_method": args.mapping_method,
        "recommendation": {
            **asdict(recommendation),
            "offsets": {str(k): v for k, v in recommendation.offsets.items()},
        },
        "read_lengths": [asdict(score) for score in scores],
    }
    (args.output_dir_path / "tis_recommendation.json").write_text(json.dumps(payload, indent=2))

    if recommendation.has_recommendation:
        print(f"{library}: read lengths {recommendation.read_lengths}, "
              f"offsets {recommendation.offsets}, {recommendation.confidence} confidence")
    else:
        print(f"{library}: no usable initiation signal, no recommendation")


if __name__ == "__main__":
    main()
