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
from lib.alignment import IntervalReader, LengthCounter


def parse_arguments(argv=None):
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
    parser.add_argument("--mapping_methods", nargs="+", default=["fiveprime", "threeprime"],
                        choices=["fiveprime", "threeprime"],
                        help="Read ends to evaluate. Both are analysed by default and the "
                             "one with the stronger initiation signal is recommended, "
                             "because which end is sharper is organism and protocol "
                             "dependent.")
    parser.add_argument("--positions_out_ORF", type=int, default=100,
                        help="Nucleotides upstream of the start codon to profile.")
    parser.add_argument("--positions_in_ORF", type=int, default=150,
                        help="Nucleotides inside the ORF to profile.")
    parser.add_argument("--length_cutoff", type=int, default=50,
                        help="Minimum ORF length when the length filter is enabled.")
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
    return parser.parse_args(argv)


def build_profiles(args):
    """Start and stop metagene profiles per read end, plus library read totals.

    The alignment is read once and reused for every read end; the annotation
    filtering and metagene mapping are redone per end, because both depend on
    which end of a read is being counted.
    """
    genome_lengths = io.parse_genome_lengths(args.genome_file_path)
    read_lengths = io.parse_read_lengths(args.read_lengths)

    reader = IntervalReader(args.alignment_file_path)
    read_intervals, total_counts = reader.output()

    profiles_by_end = {}
    for mapping_method in args.mapping_methods:
        start_codons, stop_codons = ann.retrieve_annotation_positions(
            args.annotation_file_path,
            read_intervals,
            total_counts,
            genome_lengths,
            args.filtering_methods,
            mapping_method,
            args.rpkm_threshold,
            args.neighboring_genes_distance,
            args.positions_out_ORF,
            args.positions_in_ORF,
            args.length_cutoff,
        )

        start_coverage = mg.metagene_mapping_start(
            start_codons, read_intervals, args.positions_out_ORF,
            args.positions_in_ORF, mapping_method
        )
        stop_coverage = mg.metagene_mapping_stop(
            stop_codons, read_intervals, args.positions_out_ORF,
            args.positions_in_ORF, mapping_method
        )
        start_coverage, stop_coverage = misc.equalize_dictionary_keys(
            start_coverage,
            stop_coverage,
            args.positions_out_ORF,
            args.positions_in_ORF,
            read_lengths=read_lengths,
        )

        # Sum across sequences: the offset is a property of the protocol, not of
        # a particular replicon, and pooling gives the estimate more to work with.
        window_length = args.positions_out_ORF + args.positions_in_ORF
        profiles_by_end[mapping_method] = (
            _sum_over_chromosomes(start_coverage, read_lengths, window_length),
            _sum_over_chromosomes(stop_coverage, read_lengths, window_length),
        )

    # Abundance is a property of the library, so it is counted from the reads
    # themselves rather than derived from one end's profiles.
    length_counts = LengthCounter(args.alignment_file_path, read_lengths).output()
    totals = {}
    for chromosome in length_counts:
        for read_length, count in length_counts[chromosome].items():
            totals[int(read_length)] = totals.get(int(read_length), 0) + count

    return profiles_by_end, totals


def _sum_over_chromosomes(coverage, wanted_lengths, window_length):
    """Pool contigs while retaining every explicitly requested read length."""
    pooled = {
        int(read_length): np.zeros(window_length, dtype=float)
        for read_length in wanted_lengths
    }
    for chromosome in coverage:
        for read_length, values in coverage[chromosome].items():
            key = int(read_length)
            if key not in pooled:
                continue
            values = np.asarray(values, dtype=float)
            pooled[key] += values
    return pooled


def to_dataframe(profiles, coordinates):
    import pandas as pd

    frame = pd.DataFrame({"coordinates": coordinates})
    for read_length in sorted(profiles):
        frame[str(read_length)] = profiles[read_length]
    return frame


def metagene_coordinates(positions_out_ORF, positions_in_ORF):
    """Return the distinct transcript-oriented axes for start and stop profiles."""
    start = np.arange(-positions_out_ORF, positions_in_ORF)
    stop = np.arange(-positions_in_ORF, positions_out_ORF)
    return start, stop


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------


def orfbounder_config(recommendation):
    """A config block that can be pasted straight into an ORFBounder run."""
    if not recommendation.has_recommendation:
        return "# No usable initiation signal was found; no setup is suggested."
    read_end = "5'" if recommendation.read_end == "fiveprime" else "3'"
    offsets = "\n".join(
        f"    {length}: {offset}" for length, offset in recommendation.offset_table()
    )
    return (
        "# Suggested ORFBounder settings, derived from this library.\n"
        f"readLengths: [{', '.join(str(length) for length in recommendation.read_lengths)}]\n"
        f"mappingMethod: \"{recommendation.read_end}\"\n"
        "# Offsets are measured from the "
        f"{read_end} end of the read.\n"
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


def render_report(library, best, comparisons, figures, path, include_plotly_js="integrated"):
    """The human-facing report: verdict first, then the evidence behind it."""
    recommendation = best.recommendation if best else comparisons[0].recommendation
    colour = theme.CONFIDENCE_COLORS.get(recommendation.confidence, theme.INK_MUTED)

    parts = ['<div class="verdict">']
    if recommendation.has_recommendation:
        lengths = ", ".join(str(length) for length in recommendation.read_lengths)
        end_label = "5'" if recommendation.read_end == "fiveprime" else "3'"
        parts.append(
            f"<strong>Use {end_label} mapping with read lengths {lengths}</strong>"
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

    for line in psite.describe_end_choice(best, comparisons):
        parts.append(f"<p>{line}</p>")

    for warning in recommendation.warnings:
        parts.append(f'<p class="warn">{warning}</p>')

    parts.append("<h2>Suggested configuration</h2>")
    parts.append(f"<pre><code>{orfbounder_config(recommendation)}</code></pre>")

    parts.append("<h2>5' against 3' mapping</h2>")
    parts.append(_end_comparison_table(comparisons, best))

    if recommendation.rationale:
        parts.append("<h2>How the read lengths were chosen</h2><ul>")
        parts.extend(f"<li>{line}</li>" for line in recommendation.rationale)
        parts.append("</ul>")

    for comparison in comparisons:
        end_label = "5'" if comparison.read_end == "fiveprime" else "3'"
        parts.append(f"<h2>Evidence per read length, {end_label} mapping</h2>")
        parts.append(_scores_table(comparison.scores))

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

    ends = " and ".join(c.read_end for c in comparisons)
    Path(path).write_text(
        theme.page(
            f"TIS caller advice: {library}",
            f"Both {ends} mapping were evaluated. Offsets are the distance from the "
            "mapped read end to the P-site.",
            "\n".join(parts),
        )
    )


def _end_comparison_table(comparisons, best):
    """Side by side summary, so the losing end is visible rather than discarded."""
    rows = [
        "<div class='table-wrap'><table><thead><tr>"
        "<th>Read end</th><th>Chosen</th><th>Read lengths</th><th>Offsets</th>"
        "<th>Peak vs background</th><th>Dominant frame</th><th>Periodicity</th>"
        "<th>Library covered</th><th>Confidence</th></tr></thead><tbody>"
    ]
    for comparison in comparisons:
        r = comparison.recommendation
        chosen = "yes" if best is not None and comparison.read_end == best.read_end else ""
        if r.has_recommendation:
            lengths = ", ".join(str(length) for length in r.read_lengths)
            offsets = ", ".join(f"{length}:{offset}" for length, offset in r.offset_table())
        else:
            lengths = offsets = "-"
        rows.append(
            "<tr>"
            f"<td>{comparison.read_end}</td>"
            f"<td>{chosen}</td>"
            f"<td>{lengths}</td>"
            f"<td>{offsets}</td>"
            f"<td class='num'>{r.sharpness:.1f}x</td>"
            f"<td class='num'>{r.frame_bias:.0%}</td>"
            f"<td class='num'>{r.periodicity:.0%}</td>"
            f"<td class='num'>{r.covered_fraction:.0%}</td>"
            f"<td>{r.confidence}</td>"
            "</tr>"
        )
    rows.append("</tbody></table></div>")
    return "\n".join(rows)


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

    start_coordinates, stop_coordinates = metagene_coordinates(
        args.positions_out_ORF, args.positions_in_ORF
    )
    profiles_by_end, totals = build_profiles(args)

    start_by_end = {end: start for end, (start, _) in profiles_by_end.items()}
    best, comparisons = psite.compare_read_ends(
        start_by_end, start_coordinates, totals
    )

    # Figures follow the chosen end; the other end's numbers stay in the tables.
    chosen = best.read_end if best else args.mapping_methods[0]
    start_profiles, stop_profiles = profiles_by_end[chosen]
    chosen_comparison = next(c for c in comparisons if c.read_end == chosen)
    scores = chosen_comparison.scores
    recommendation = chosen_comparison.recommendation

    read_lengths = sorted(start_profiles)
    df_start = to_dataframe(start_profiles, start_coordinates)
    df_stop = to_dataframe(stop_profiles, stop_coordinates)
    significant_offsets = {s.read_length: (s.offset if s.usable else None) for s in scores}
    end_label = "5'" if chosen == "fiveprime" else "3'"

    figures = [
        ("Read length against position", plotting.plot_metagene_heatmap(
            df_start, df_stop, read_lengths, "Metagene profile",
            f"{end_label} mapping, enrichment over each read length's own background",
            significant_offsets, read_end=chosen)),
        ("Library composition", plotting.plot_read_length_distribution(
            scores, "Read length distribution",
            f"Which read lengths carry a usable initiation signal under {end_label} mapping")),
        ("Reading frame", plotting.plot_frame_composition(
            scores, "Reading frame composition",
            "Share of P-sites per frame, after offset correction")),
    ]
    if recommendation.has_recommendation:
        pooled = psite.pool_profiles(
            start_profiles, recommendation.offsets, recommendation.read_lengths,
            psite.READ_ENDS[chosen]
        )
        figures.append(("Pooled, offset corrected", plotting.plot_pooled_profile(
            pooled, start_coordinates, "Pooled P-site profile",
            f"{end_label} mapping, read lengths "
            f"{', '.join(str(length) for length in recommendation.read_lengths)}")))

    render_report(library, best, comparisons, figures,
                  args.output_dir_path / "tis_recommendation.html", args.include_plotly_js)
    scores_to_tsv(scores, args.output_dir_path / "read_length_evidence.tsv")

    payload = {
        "library": library,
        "evaluated_read_ends": list(profiles_by_end),
        "chosen_read_end": chosen if best else None,
        "recommendation": {
            **asdict(recommendation),
            "offsets": {str(k): v for k, v in recommendation.offsets.items()},
        },
        "read_ends": {
            comparison.read_end: {
                "quality": comparison.quality,
                "recommendation": {
                    **asdict(comparison.recommendation),
                    "offsets": {str(k): v for k, v in comparison.recommendation.offsets.items()},
                },
                "read_lengths": [asdict(score) for score in comparison.scores],
            }
            for comparison in comparisons
        },
    }
    (args.output_dir_path / "tis_recommendation.json").write_text(json.dumps(payload, indent=2))

    if recommendation.has_recommendation:
        print(f"{library}: {chosen} mapping, read lengths {recommendation.read_lengths}, "
              f"offsets {recommendation.offsets}, {recommendation.confidence} confidence")
    else:
        print(f"{library}: no usable initiation signal on either read end, no recommendation")


if __name__ == "__main__":
    main()
