#!/usr/bin/env python
"""
Recommend read lengths and site offsets for initiation or termination peaks.

Given an alignment and an annotation, this works out which read lengths carry a
usable boundary signal, how far each of them sits from the occupied site, and whether
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
import html
from collections import defaultdict
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
from lib.orfbounder import export_orfbounder_inputs, recommendation_inputs


DEEPRIBO_DEFAULT_A_SITE_OFFSET = 12
DEEPRIBO_MIN_AGREEMENT = 0.80
DEEPRIBO_MIN_READ_SUPPORT = 0.50
DEEPRIBO_SINGLE_LENGTH_READ_SUPPORT = 0.80


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Recommend read lengths and P-site (TIS) or A-site (TTS) offsets.",
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
    parser.add_argument("--library_type", choices=["auto", "RIBO", "TIS", "TTS"],
                        default="auto", help="TTS evaluates stop peaks and A-site offsets; "
                        "RIBO/TIS evaluate start peaks and P-site offsets. Auto uses "
                        "the BAM filename's method prefix, defaulting to RIBO.")
    parser.add_argument("--mapping_methods", nargs="+", default=["fiveprime", "threeprime"],
                        choices=["fiveprime", "threeprime"],
                        help="Read ends to evaluate. Both are analysed by default and the "
                             "one with the more consistent boundary offsets is recommended, "
                             "because which end is sharper is organism and protocol "
                             "dependent.")
    parser.add_argument("--positions_out_ORF", type=int, default=100,
                        help="Nucleotides outside each ORF boundary to profile.")
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
    parser.add_argument("--deepribo_asite_offset", type=int,
                        default=DEEPRIBO_DEFAULT_A_SITE_OFFSET,
                        help="Current DeepRibo 3'-to-A-site offset, shown for comparison only.")
    parser.add_argument("--min_relative_enrichment", type=float,
                        default=psite.DEFAULT_MIN_RELATIVE_ENRICHMENT,
                        help="Minimum fraction of the coverage-supported reference's "
                        "enrichment retained by each length and the pooled profile (0 to 1). "
                        "Lower values favor broader coverage; higher values favor peak quality.")
    args = parser.parse_args(argv)
    if not np.isfinite(args.min_relative_enrichment) or not 0 <= args.min_relative_enrichment <= 1:
        parser.error("--min_relative_enrichment must be a finite number between 0 and 1")
    return args


def build_profiles(args):
    """Start and stop metagene profiles per read end, plus library read totals.

    Alignment and annotation filtering are shared by both read ends, so their
    offset consistency is compared on exactly the same retained genes.
    """
    genome_lengths = io.parse_genome_lengths(args.genome_file_path)
    read_lengths = io.parse_read_lengths(args.read_lengths)

    reader = IntervalReader(args.alignment_file_path)
    read_intervals, total_counts = reader.output()

    # Initiation 5' ends lie before the CDS; termination 3' ends lie after it.
    # Gene abundance belongs to overlapping footprints, not a selected end.
    start_codons, stop_codons = ann.retrieve_annotation_positions(
        args.annotation_file_path,
        read_intervals,
        total_counts,
        genome_lengths,
        args.filtering_methods,
        "global",
        args.rpkm_threshold,
        args.neighboring_genes_distance,
        args.positions_out_ORF,
        args.positions_in_ORF,
        args.length_cutoff,
    )

    profiles_by_end = {}
    for mapping_method in args.mapping_methods:
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

    return profiles_by_end, totals, sum(total_counts.values())


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


def library_type_from_name(library, requested="auto"):
    """Resolve workflow sample methods and standalone BAM overrides."""
    if requested != "auto":
        return requested
    method = library.split("-", 1)[0].upper()
    return method if method in {"RIBO", "TIS", "TTS"} else "RIBO"


def site_offsets(offsets, read_end, site):
    """Distances to the first bases of adjacent P- and A-site codons."""
    direction = psite.READ_ENDS[read_end].sign
    if site == "A":
        a_offsets = dict(offsets)
        p_offsets = {length: offset - direction * 3 for length, offset in offsets.items()}
    else:
        p_offsets = dict(offsets)
        a_offsets = {length: offset + direction * 3 for length, offset in offsets.items()}
    return p_offsets, a_offsets


def recommendation_payload(recommendation):
    p_offsets, a_offsets = site_offsets(
        recommendation.offsets, recommendation.read_end, recommendation.site
    )
    return {
        **asdict(recommendation),
        "offsets": {str(k): v for k, v in recommendation.offsets.items()},
        "p_site_offsets": {str(k): v for k, v in p_offsets.items()},
        "a_site_offsets": {str(k): v for k, v in a_offsets.items()},
    }


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------


def orfbounder_config(recommendation, library):
    """Preview the two actual ORFBounder JSON inputs and their mapping method."""
    if not recommendation.has_recommendation:
        return "No supported calibration; input JSON files are not exported."
    lengths, offsets = recommendation_inputs(
        library.split("_", 1)[0], recommendation_payload(recommendation), recommendation.read_end
    )
    return (
        f"mapping_method: {recommendation.read_end}\n\n"
        f"read_lengths.json:\n{json.dumps(lengths, indent=2)}\n\n"
        f"offsets.json:\n{json.dumps(offsets, indent=2)}\n"
    )


def deepribo_asite_advice(library, comparisons, total_reads, current_offset,
                         library_type=None):
    """Suggest a DeepRibo scalar only when usable 3' peaks agree across reads.

    DeepRibo currently uses one 3'-to-A-site offset for the RIBO reads it accepts.
    A site's first base is three nucleotides downstream of the P site's first
    base, so the candidate from a 3'-to-P offset is offset - 3. The TIS caller's
    greedy read-length subset is intentionally not used here: it can conceal
    conflicting estimates at other abundant lengths that DeepRibo still reads.
    """
    advice = {
        "current_offset": current_offset,
        "suggested_offset": None,
        "applied": False,
        "per_length_offsets": {},
        "supporting_read_lengths": [],
        "agreement_fraction": 0.0,
        "read_support_fraction": 0.0,
        "reason": "",
    }
    if (library_type is not None and library_type != "RIBO") or (
        library_type is None and not library.startswith("RIBO-")
    ):
        advice["reason"] = "DeepRibo uses RIBO libraries; this is not a RIBO library."
        return advice

    threeprime = next((c for c in comparisons if c.read_end == "threeprime"), None)
    if threeprime is None:
        advice["reason"] = "3' read ends were not evaluated."
        return advice
    if threeprime.recommendation.confidence not in {"medium", "high"}:
        advice["reason"] = "The 3' initiation signal is not strong enough for a DeepRibo suggestion."
        return advice
    fractions = threeprime.recommendation.frame_fractions
    if max(fractions) >= 0.5 and fractions[0] < max(fractions):
        advice["reason"] = (
            "The dominant 3' P-site reading frame is not frame 0; the estimated "
            "offsets may be shifted, so no DeepRibo value is suggested."
        )
        return advice

    by_offset = defaultdict(list)
    for score in threeprime.scores:
        if not score.usable or score.offset is None or score.total_reads <= 0:
            continue
        candidate = score.offset - 3
        if not 0 <= candidate < score.read_length or candidate > 30:
            continue
        advice["per_length_offsets"][str(score.read_length)] = candidate
        by_offset[candidate].append(score)

    if not by_offset or total_reads <= 0:
        advice["reason"] = "No usable 3' P-site estimates yield a valid A-site offset."
        return advice

    weights = {
        offset: sum(score.total_reads for score in scores)
        for offset, scores in by_offset.items()
    }
    max_weight = max(weights.values())
    leaders = [offset for offset, weight in weights.items() if weight == max_weight]
    if len(leaders) != 1:
        advice["reason"] = "Usable read lengths have equally supported but different A-site offsets."
        return advice

    candidate = leaders[0]
    supporting_scores = by_offset[candidate]
    usable_reads = sum(weights.values())
    agreement = max_weight / usable_reads
    read_support = max_weight / total_reads
    advice["supporting_read_lengths"] = sorted(s.read_length for s in supporting_scores)
    advice["agreement_fraction"] = round(agreement, 4)
    advice["read_support_fraction"] = round(read_support, 4)

    if agreement < DEEPRIBO_MIN_AGREEMENT:
        advice["reason"] = (
            "Usable 3' read lengths disagree on a single A-site offset "
            f"({agreement:.0%} agreement; at least {DEEPRIBO_MIN_AGREEMENT:.0%} required)."
        )
    elif read_support < DEEPRIBO_MIN_READ_SUPPORT:
        advice["reason"] = (
            "The best A-site offset represents too few reads accepted by the advisor "
            f"({read_support:.0%}; at least {DEEPRIBO_MIN_READ_SUPPORT:.0%} required)."
        )
    elif (
        len(supporting_scores) == 1
        and read_support < DEEPRIBO_SINGLE_LENGTH_READ_SUPPORT
    ):
        advice["reason"] = (
            "Only one read length supports the A-site offset, without dominating "
            f"the library ({read_support:.0%}; at least "
            f"{DEEPRIBO_SINGLE_LENGTH_READ_SUPPORT:.0%} required)."
        )
    else:
        advice["suggested_offset"] = candidate
        advice["reason"] = (
            f"{agreement:.0%} of usable 3' reads and {read_support:.0%} of all "
            "reads accepted by the advisor support this offset."
        )
    return advice


def scores_to_tsv(scores, path, read_end="fiveprime", anchor="start", site="P"):
    header = [
        "read_length", "total_reads", "abundance", "offset", "peak_height",
        "background", "sharpness", "z_score", "frame_0", "frame_1", "frame_2",
        "frame_bias", "periodicity", "usable", "reasons",
        "read_end", "anchor", "site", "p_site_offset", "a_site_offset",
        "background_reference", "background_positions",
        "frame_reads",
    ]
    with open(path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for score in scores:
            p_offsets, a_offsets = site_offsets(
                {score.read_length: score.offset} if score.offset is not None else {},
                read_end, site,
            )
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
                read_end, anchor, site,
                str(p_offsets.get(score.read_length, "")),
                str(a_offsets.get(score.read_length, "")),
                f"{score.background_reference:.6f}", str(score.background_positions),
                "" if score.frame_reads is None else f"{score.frame_reads:.0f}",
            ]) + "\n")


def render_report(library, best, comparisons, figures, asite_advice, path,
                  include_plotly_js="integrated", library_type="RIBO"):
    """The human-facing report: verdict first, then the evidence behind it."""
    recommendation = best.recommendation if best else comparisons[0].recommendation
    anchor = recommendation.anchor
    site = recommendation.site
    signal = "termination" if anchor == "stop" else "initiation"
    advice_label = "TTS peak advice" if anchor == "stop" else "TIS caller advice"
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
            f"<p>The pooled {signal} peak is {recommendation.sharpness:.1f} times the "
            f"upstream background, built from {recommendation.covered_fraction:.1%} of the "
            f"library. The dominant reading frame holds {recommendation.frame_bias:.0%} of "
            f"{site}-sites.</p>"
        )
    else:
        parts.append(
            f"<strong>No {'TTS peak setup' if anchor == 'stop' else 'TIS caller setup'} "
            "is recommended for this library</strong>"
            f'<span class="pill" style="background:{colour}">no recommendation</span>'
        )
    parts.append("</div>")

    for line in psite.describe_end_choice(best, comparisons, anchor=anchor):
        parts.append(f"<p>{line}</p>")

    for warning in recommendation.warnings:
        parts.append(f'<p class="warn">{warning}</p>')

    parts.append("<h2>Site offsets</h2>")
    parts.append(
        f"<p>This {library_type} library is scored at the {anchor} codon. "
        f"Reported offsets measure the distance to the first nucleotide of the {site}-site codon "
        "from the inclusive mapped read end. From a 5' end, move downstream; "
        "from a 3' end, move upstream, in transcript direction on either strand.</p>"
    )
    if anchor == "stop":
        parts.append(
            "<p>Apidaecin termination complexes place the stop codon in the A-site. "
            "The P-site is one codon (3 nt) upstream: a 5'-to-P-site distance is "
            "the A-site offset minus 3; a 3'-to-P-site distance is the A-site offset plus 3. "
            "These are termination-peak settings.</p>"
            "<p>Interpret termination enrichment with the protocol in mind: apidaecin can "
            "also produce upstream ribosome queues and stop-codon readthrough. "
            "For disome footprints, include their lengths in the configured range and "
            "interpret the estimate as the leading terminating ribosome's site.</p>"
            '<p>Scientific basis: <a href="https://www.nature.com/articles/s41467-025-58329-w">'
            'Froschauer et al. (2025)</a> and '
            '<a href="https://elifesciences.org/articles/62655">Mangano et al. (2020)</a>.</p>'
        )
    parts.append(_site_offsets_table(recommendation))
    parts.append("<h2>ORFBounder JSON inputs</h2>")
    parts.append(
        "<p>The separate read-length and offset JSON files use this library's "
        "supported calibration for each mapped end. ORFBounder uses one "
        "<code>mapping_method</code> for every assay in a run. Exported 5' offsets "
        "are negative and 3' offsets positive, aligning to the first nucleotide "
        f"of the calibrated {site}-site codon. These files prepare inputs; they "
        "do not run ORFBounder.</p>"
        '<p><a href="orfbounder/manifest.json">Export availability and calibration metadata</a></p>'
    )
    for comparison in comparisons:
        proposed = comparison.recommendation
        if not proposed.has_recommendation:
            continue
        end = comparison.read_end
        parts.append(f"<h3>{end}" + (" (recommended end)" if best and end == best.read_end else "") + "</h3>")
        parts.append(
            f'<p><a href="orfbounder/{end}/read_lengths.json">read_lengths.json</a> · '
            f'<a href="orfbounder/{end}/offsets.json">offsets.json</a></p>'
        )
        parts.append(f"<pre><code>{html.escape(orfbounder_config(proposed, library))}</code></pre>")

    if anchor == "start":
        parts.append("<h2>DeepRibo A-site offset advice</h2>")
        suggested = asite_advice["suggested_offset"]
        current = asite_advice["current_offset"]
        if suggested is None:
            parts.append("<p>No single DeepRibo offset is suggested for this library.</p>")
        elif suggested == current:
            parts.append(f"<p>The suggested DeepRibo A-site offset is {suggested} nt, "
                         "matching the current configuration.</p>")
        else:
            parts.append(f"<p>The suggested DeepRibo A-site offset is {suggested} nt; "
                         f"the current configuration is {current} nt.</p>")
            parts.append("<p>In the existing <code>predictionSettings</code> block, "
                         "change only this line if the advice is appropriate:</p>")
            parts.append(f"<pre><code>  deepriboASiteOffset: {suggested}</code></pre>")
        parts.append(f"<p>{asite_advice['reason']}</p>")
        if asite_advice["per_length_offsets"]:
            values = ", ".join(
                f"{length} nt → {offset} nt"
                for length, offset in sorted(
                    asite_advice["per_length_offsets"].items(),
                    key=lambda item: int(item[0]),
                )
            )
            parts.append(f"<p>Usable 3' read-length candidates: {values}.</p>")
        parts.append(
            "<p>This is advice only: HRIBO does not change DeepRibo automatically. "
            "The estimate assumes the A-site is one codon downstream of the P-site. "
            "DeepRibo uses one global offset, so compare advice across RIBO "
            "libraries before changing the setting and rerunning predictions.</p>"
        )

    parts.append("<h2>5' against 3' mapping</h2>")
    parts.append(_end_comparison_table(comparisons, best, site=site))

    if recommendation.rationale:
        parts.append("<h2>How the read lengths were chosen</h2><ul>")
        parts.extend(f"<li>{line}</li>" for line in recommendation.rationale)
        parts.append("</ul>")
    parts.append(
        "<p>Selection uses a reference ranked by read share times log(1 + enrichment), "
        "then adds the largest available share of reads while each length and the pooled "
        "profile meet the reported quality floor and the noise/background checks. "
        "An addition can reduce enrichment while retaining a strong pooled peak. "
        "The configurable relative floor is an analysis preference, not a validated "
        "biological cutoff; the search adds one feasible length at a time. "
        "Raw peak depth, enrichment, and coding-body reading frames are separate measurements. "
        "The evidence tables and scored heatmap panel use the same upstream background "
        "relative to the calibrated site; the other panel is a diagnostic normalized "
        "by its whole-window median. Hover over the heatmap for exact enrichment and read counts.</p>"
    )
    parts.append(
        "<p>The 3-nt FFT score describes the coding-body metagene away from the boundary peak. "
        "It measures the share of non-DC FFT power in the bin nearest a three-nucleotide period. "
        "It is not a percentage of correctly assigned reads or a statistical significance test. "
        "A low score alone does not invalidate bacterial start- or stop-codon offsets; "
        "window length, coverage shape, and nuclease bias affect this diagnostic.</p>"
    )

    for comparison in comparisons:
        end_label = "5'" if comparison.read_end == "fiveprime" else "3'"
        parts.append(f"<h2>Evidence per read length, {end_label} mapping</h2>")
        parts.append(_scores_table(comparison.scores, site=site,
                                   selected_lengths=comparison.recommendation.read_lengths))

    parts.append("""<script>
(() => {
  document.querySelectorAll('.adviser-table-wrap').forEach(wrapper => {
    const hint = document.createElement('p');
    hint.className = 'adviser-scroll-hint';
    hint.textContent = 'Scroll sideways to see all columns →';
    hint.hidden = true;
    wrapper.before(hint);
    const update = () => {
      const overflowing = wrapper.scrollWidth > wrapper.clientWidth + 1;
      hint.hidden = !overflowing;
      if (overflowing) {
        wrapper.tabIndex = 0;
        wrapper.setAttribute('role', 'region');
        wrapper.setAttribute('aria-label', 'Scrollable table; use left and right arrow keys');
      } else {
        wrapper.removeAttribute('tabindex');
        wrapper.removeAttribute('role');
        wrapper.removeAttribute('aria-label');
      }
    };
    update();
    window.addEventListener('resize', update);
    if (typeof ResizeObserver !== 'undefined') {
      const observer = new ResizeObserver(update);
      observer.observe(wrapper);
      observer.observe(wrapper.querySelector('table'));
    }
  });
})();
</script>""")

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
            f"{advice_label}: {library}",
            f"Evaluated read ends: {ends}. {signal.capitalize()} peaks at {anchor} codons; "
            f"offsets run from the mapped read end to the {site}-site.",
            "\n".join(parts),
        )
    )


def _site_offsets_table(recommendation):
    p_offsets, a_offsets = site_offsets(
        recommendation.offsets, recommendation.read_end, recommendation.site
    )
    if recommendation.site == "P":
        measured_site, derived_site = "P", "A"
        measured_offsets, derived_offsets = p_offsets, a_offsets
    else:
        measured_site, derived_site = "A", "P"
        measured_offsets, derived_offsets = a_offsets, p_offsets
    rows = ["<div class='table-wrap adviser-table-wrap'><table class='adviser-table adviser-offsets'>"
            "<thead><tr><th class='num'>Read length</th>"
            f"<th class='num'>{measured_site}-site offset (nt)</th>"
            f"<th class='num'>Derived {derived_site}-site offset (nt)</th>"
            "</tr></thead><tbody>"]
    rows.extend(
        f"<tr><td class='num'>{length}</td><td class='num'>{measured_offsets[length]}</td>"
        f"<td class='num'>{derived_offsets[length]}</td></tr>"
        for length in sorted(recommendation.read_lengths)
    )
    rows.append("</tbody></table></div>")
    return "\n".join(rows)


def _end_comparison_table(comparisons, best, site="P"):
    """Side by side summary, so the losing end is visible rather than discarded."""
    rows = [
        "<div class='table-wrap adviser-table-wrap'><table class='adviser-table'><thead><tr>"
        f"<th>Read end</th><th>Chosen</th><th>Read lengths</th><th>{site}-site offsets (nt)</th>"
        "<th>Peak vs background</th><th>Dominant frame</th><th>3-nt FFT score</th>"
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


def _scores_table(scores, site="P", selected_lengths=()):
    selected = set(selected_lengths)
    rows = [
        "<div class='table-wrap adviser-table-wrap'><table class='adviser-table adviser-evidence'><thead><tr>"
        f"<th class='num'>Read length</th><th class='num'>Share of reads</th>"
        f"<th class='num'>{site}-site offset (nt)</th>"
        "<th class='num'>Peak read ends</th><th class='num'>Background read ends/bin</th>"
        "<th class='num'>Peak vs background</th><th class='num'>Body read ends</th>"
        "<th class='num'>Frame 0</th><th class='num'>Dominant frame</th>"
        "<th class='num'>3-nt FFT score</th>"
        "<th>Usable</th><th>Selected</th><th class='notes'>Notes</th></tr></thead><tbody>"
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
            f"<td class='num'>{score.peak_height:.0f}</td>"
            f"<td class='num'>{score.background_reference:.2f}</td>"
            f"<td class='num'>{score.sharpness:.1f}x</td>"
            f"<td class='num'>{'-' if score.frame_reads is None else f'{score.frame_reads:.0f}'}</td>"
            f"<td class='num'>{score.frame_fractions[0]:.0%}</td>"
            f"<td class='num'>{frame}</td>"
            f"<td class='num'>{score.periodicity:.0%}</td>"
            f"<td>{'yes' if score.usable else 'no'}</td>"
            f"<td>{'yes' if score.read_length in selected else 'no'}</td>"
            f"<td class='notes'>{'; '.join(score.reasons)}</td>"
            "</tr>"
        )
    rows.append("</tbody></table></div>")
    return "\n".join(rows)


def main():
    args = parse_arguments()
    args.output_dir_path.mkdir(parents=True, exist_ok=True)
    library = args.alignment_file_path.stem
    library_type = library_type_from_name(library, args.library_type)
    anchor = "stop" if library_type == "TTS" else "start"
    site = "A" if anchor == "stop" else "P"
    signal = "termination" if anchor == "stop" else "initiation"

    start_coordinates, stop_coordinates = metagene_coordinates(
        args.positions_out_ORF, args.positions_in_ORF
    )
    # The legacy stop axis places the three stop bases at -3, -2, -1.
    # Offset distances must use the FIRST stop base, just like start distances.
    stop_coordinates = stop_coordinates + 3
    profiles_by_end, totals, total_reads = build_profiles(args)

    boundary_by_end = {
        end: stop if anchor == "stop" else start
        for end, (start, stop) in profiles_by_end.items()
    }
    coordinates = stop_coordinates if anchor == "stop" else start_coordinates
    best, comparisons = psite.compare_read_ends(
        boundary_by_end, coordinates, totals, anchor=anchor,
        min_relative_enrichment=args.min_relative_enrichment,
    )
    asite_advice = deepribo_asite_advice(
        library, comparisons, total_reads, args.deepribo_asite_offset, library_type
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
            significant_offsets, read_end=chosen, anchor=anchor, site=site,
            background_references={s.read_length: s.background_reference for s in scores})),
        ("Library composition", plotting.plot_read_length_distribution(
            scores, "Read length distribution",
            f"Which read lengths carry a usable {signal} signal under {end_label} mapping",
            anchor=anchor)),
        ("Reading frame", plotting.plot_frame_composition(
            scores, "Reading frame composition",
            f"Share of {site}-sites per frame inside the ORF, after offset correction",
            site=site)),
    ]
    if recommendation.has_recommendation:
        pooled = psite.pool_profiles(
            boundary_by_end[chosen], recommendation.offsets, recommendation.read_lengths,
            psite.READ_ENDS[chosen]
        )
        figures.append(("Pooled, offset corrected", plotting.plot_pooled_profile(
            pooled, coordinates, f"Pooled {site}-site profile",
            f"{end_label} mapping, read lengths "
            f"{', '.join(str(length) for length in recommendation.read_lengths)}",
            anchor=anchor, site=site)))

    render_report(library, best, comparisons, figures, asite_advice,
                  args.output_dir_path / "tis_recommendation.html", args.include_plotly_js,
                  library_type=library_type)
    scores_to_tsv(scores, args.output_dir_path / "read_length_evidence.tsv",
                  chosen, anchor, site)

    payload = {
        "library": library,
        "library_type": library_type,
        "anchor": anchor,
        "site": site,
        "offset_reference": "first nucleotide of the site codon, from the inclusive mapped read end",
        "evaluated_read_ends": list(profiles_by_end),
        "chosen_read_end": chosen if best else None,
        "deepribo_a_site": asite_advice,
        "recommendation": recommendation_payload(recommendation),
        "read_ends": {
            comparison.read_end: {
                "quality": comparison.quality,
                "recommendation": recommendation_payload(comparison.recommendation),
                "read_lengths": [asdict(score) for score in comparison.scores],
            }
            for comparison in comparisons
        },
    }
    (args.output_dir_path / "tis_recommendation.json").write_text(json.dumps(payload, indent=2))
    export_orfbounder_inputs([payload], args.output_dir_path / "orfbounder")

    if recommendation.has_recommendation:
        print(f"{library}: {chosen} mapping, read lengths {recommendation.read_lengths}, "
              f"{site}-site offsets {recommendation.offsets}, {recommendation.confidence} confidence")
    else:
        print(f"{library}: no usable {signal} signal on the evaluated read ends, no recommendation")


if __name__ == "__main__":
    main()
