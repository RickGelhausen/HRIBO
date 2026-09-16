#!/usr/bin/env python3
"""Summarize condition-level detection and differential results across contrasts.

The two questions are deliberately kept separate. Detection is a descriptive,
replicate-aware threshold on within-assay counts per million (CPM); differential
calls come from deltaTE's RNA, RIBO, and TE models. xTail TE values, and optional
RiboRex values when that analysis is enabled, remain supplementary estimates
and are never combined into a new p-value.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

import xlsxwriter

from lib.diffex_browser import export_tracks


ASSAYS = ("RNA", "RIBO")
EFFECT_ASSAYS = ("RNA", "RIBO", "TE")
CONDITION_FIELDS = (
    "feature_id", "locus_tag", "genome", "start", "end", "strand",
    "feature_type", "name", "condition", "assay", "state", "replicates",
    "passing_replicates", "mean_count", "mean_cpm", "sample_counts",
    "sample_cpms",
)
CONTRAST_FIELDS = (
    "feature_id", "locus_tag", "genome", "start", "end", "strand",
    "feature_type", "name", "contrast", "assay", "state", "log2fc", "padj",
    "method", "xtail_te_log2fc", "xtail_te_padj", "riborex_te_log2fc",
    "riborex_te_padj",
)
TRACK_FIELDS = (
    "path", "kind", "condition", "contrast", "assay", "direction",
    "feature_count",
)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--xtail", type=Path, required=True)
    parser.add_argument(
        "--riborex", type=Path,
        help="Optional pooled RiboRex results; omit when RiboRex was not run",
    )
    parser.add_argument("--deltate", type=Path, required=True)
    parser.add_argument("--contrasts", nargs="+", required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--min_cpm", type=float, default=1.0)
    parser.add_argument("--min_count", type=int, default=10)
    parser.add_argument("--min_replicates", type=int, default=2)
    parser.add_argument("--padj_cutoff", type=float, default=0.05)
    parser.add_argument("--log2fc_cutoff", type=float, default=1.0)
    args = parser.parse_args(argv)
    if not math.isfinite(args.min_cpm) or args.min_cpm < 0:
        parser.error("--min_cpm must be a finite non-negative number")
    if args.min_count < 1 or args.min_replicates < 1:
        parser.error("--min_count and --min_replicates must be positive integers")
    if not 0 < args.padj_cutoff < 1 or not math.isfinite(args.padj_cutoff):
        parser.error("--padj_cutoff must be between 0 and 1")
    if not math.isfinite(args.log2fc_cutoff) or args.log2fc_cutoff < 0:
        parser.error("--log2fc_cutoff must be a finite non-negative number")
    for contrast in args.contrasts:
        if len(contrast.split("-")) != 2 or not all(contrast.split("-")):
            parser.error(f"Invalid contrast {contrast!r}; expected LEFT-RIGHT")
    return args


def parse_sample(label):
    parts = label.split("-")
    if len(parts) != 3 or parts[0] not in ASSAYS or not parts[1] or not parts[2]:
        raise ValueError(f"Unexpected count-matrix sample column {label!r}")
    return parts[0], parts[1], parts[2]


def _number(value, context, *, nonnegative=False):
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid number for {context}: {value!r}") from error
    if not math.isfinite(number) or (nonnegative and number < 0):
        raise ValueError(f"Invalid number for {context}: {value!r}")
    return number


def read_counts(path):
    """Return feature-by-sample counts, sample labels, and assigned-count totals."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or reader.fieldnames[0] != "Identifier":
            raise ValueError(f"{path} must begin with an Identifier column")
        all_labels = reader.fieldnames[1:]
        if len(all_labels) != len(set(all_labels)):
            raise ValueError(f"{path} needs distinct RNA/RIBO sample columns")
        # featureCounts receives every mapped library, including TIS/TTS and
        # their RNA controls. Differential tools use only matched RNA/RIBO.
        labels = [label for label in all_labels if label.split("-", 1)[0] in ASSAYS]
        if not labels:
            raise ValueError(f"{path} has no RNA/RIBO sample columns")
        for label in labels:
            parse_sample(label)
        counts = {}
        totals = dict.fromkeys(labels, 0.0)
        for row in reader:
            feature_id = row["Identifier"]
            if not feature_id or feature_id in counts:
                raise ValueError(f"Missing or duplicate feature ID in {path}: {feature_id!r}")
            values = {
                label: _number(row[label], f"{feature_id}/{label}", nonnegative=True)
                for label in labels
            }
            counts[feature_id] = values
            for label, value in values.items():
                totals[label] += value
    return counts, labels, totals


def read_annotation(path):
    """Map feature-count coordinate IDs to GFF3 coordinates and display names."""
    coordinates = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}:{line_no}: expected at least nine GFF columns")
            seqid, _, feature_type, start, end, _, strand, _, attr_text = fields[:9]
            attrs = {}
            for item in attr_text.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = unquote(value)
                elif len(item.strip().split(None, 1)) == 2:
                    key, value = item.strip().split(None, 1)
                    attrs[key] = unquote(value.strip().strip('"'))
            feature_id = f"{seqid}:{start}-{end}:{strand}"
            record = {
                "seqid": seqid, "start": int(start), "end": int(end),
                "strand": strand, "type": feature_type,
                "locus_tag": attrs.get("locus_tag") or "",
                "gene_id": attrs.get("gene_id") or "",
                "name": (attrs.get("Name") or attrs.get("gene_name")
                         or attrs.get("locus_tag") or attrs.get("gene_id")
                         or attrs.get("ID") or feature_id),
                "aliases": list(dict.fromkeys(
                    attrs[key] for key in (
                        "Name", "gene_name", "locus_tag", "old_locus_tag",
                        "gene_id", "gene", "ID",
                    ) if attrs.get(key)
                )),
            }
            # The unambiguous annotation can contain a gene and a counted
            # feature at identical coordinates. Prefer the feature-level row.
            prior = coordinates.get(feature_id)
            replace_prior = prior is None or (
                prior["type"] in {"gene", "pseudogene", "exon"}
                and feature_type not in {"gene", "pseudogene", "exon"}
            )
            if replace_prior:
                if prior is not None:
                    if not record["locus_tag"]:
                        record["locus_tag"] = prior["locus_tag"]
                    if not record["gene_id"]:
                        record["gene_id"] = prior["gene_id"]
                    record["aliases"] = list(dict.fromkeys(
                        [*record["aliases"], *prior["aliases"]]
                    ))
                coordinates[feature_id] = record
            elif prior is not None:
                if not prior["locus_tag"] and record["locus_tag"]:
                    prior["locus_tag"] = record["locus_tag"]
                if not prior["gene_id"] and record["gene_id"]:
                    prior["gene_id"] = record["gene_id"]
                prior["aliases"] = list(dict.fromkeys(
                    [*prior["aliases"], *record["aliases"]]
                ))
    return coordinates


def read_effects(path, method):
    """Index pooled method results by (feature ID, unprefixed contrast)."""
    indexed = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"gene_id", "contrast"} <= set(reader.fieldnames):
            raise ValueError(f"{path} lacks gene_id or contrast")
        for row in reader:
            raw_contrast = row["contrast"]
            prefix = f"{method}_"
            if not raw_contrast.startswith(prefix):
                raise ValueError(f"Unexpected {method} contrast in {path}: {raw_contrast!r}")
            key = (row["gene_id"], raw_contrast[len(prefix):])
            if key in indexed:
                raise ValueError(f"Duplicate {method} result for {key!r}")
            indexed[key] = row
    return indexed


def optional_float(value):
    if value is None or str(value).strip() in ("", "NA", "NaN", "nan", "None"):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def effect_call(result, assay, padj_cutoff, log2fc_cutoff):
    """Classify a deltaTE component; missing/NA is not non-significance."""
    fc = optional_float(result.get(f"{assay}_log2FC")) if result else None
    padj = optional_float(result.get(f"{assay}_pvalue_adjusted")) if result else None
    if fc is None or padj is None:
        state = "not_tested"
    elif padj <= padj_cutoff and fc >= log2fc_cutoff and fc > 0:
        state = "up"
    elif padj <= padj_cutoff and fc <= -log2fc_cutoff and fc < 0:
        state = "down"
    else:
        state = "not_significant"
    return {"state": state, "log2fc": fc, "padj": padj}


def detection_call(values, labels, totals, min_cpm, min_count, min_replicates):
    """Classify counts in one condition/assay without claiming biological absence."""
    samples = []
    for label in labels:
        total = totals[label]
        count = values[label]
        cpm = count * 1_000_000 / total if total > 0 else None
        samples.append({"sample": label, "count": count, "cpm": cpm})
    usable = [sample for sample in samples if sample["cpm"] is not None]
    passing = sum(
        sample["count"] >= min_count and sample["cpm"] >= min_cpm
        for sample in usable
    )
    if len(usable) < min_replicates:
        state = "uncertain"
    elif passing >= min_replicates:
        state = "detected"
    elif passing == 0:
        state = "not_detected"
    else:
        state = "uncertain"
    mean_count = sum(sample["count"] for sample in usable) / len(usable) if usable else None
    mean_cpm = sum(sample["cpm"] for sample in usable) / len(usable) if usable else None
    return {
        "state": state, "replicates": len(usable), "passing_replicates": passing,
        "mean_count": mean_count, "mean_cpm": mean_cpm,
        "normalized_count": mean_cpm,
        "replicate_fraction": passing / len(usable) if usable else None,
        "samples": samples,
    }


def _format(value):
    return "" if value is None else format(value, ".8g") if isinstance(value, float) else str(value)


def _write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _metadata(feature_id, coordinates):
    item = coordinates.get(feature_id, {})
    return {
        "feature_id": feature_id,
        "locus_tag": item.get("locus_tag") or item.get("gene_id", ""),
        "genome": item.get("seqid", ""), "start": item.get("start", ""),
        "end": item.get("end", ""), "strand": item.get("strand", ""),
        "feature_type": item.get("type", ""),
        "name": item.get("name", feature_id),
    }


def build_summary(counts, labels, totals, coordinates, contrasts, effects, args):
    """Return export rows and two complete, tidy TSV tables."""
    sample_groups = defaultdict(list)
    for label in labels:
        assay, condition, _ = parse_sample(label)
        sample_groups[(condition, assay)].append(label)
    ordered_conditions = []
    for contrast in contrasts:
        left, right = contrast.split("-")
        for condition in (right, left):
            if condition not in ordered_conditions:
                ordered_conditions.append(condition)
    for _, condition, _ in map(parse_sample, labels):
        if condition not in ordered_conditions:
            ordered_conditions.append(condition)

    rows, condition_table, contrast_table = [], [], []
    for feature_id, values in counts.items():
        metadata = _metadata(feature_id, coordinates)
        row = {"feature_id": feature_id, "conditions": {}, "contrasts": {}}
        for condition in ordered_conditions:
            row["conditions"][condition] = {}
            for assay in ASSAYS:
                call = detection_call(
                    values, sample_groups[(condition, assay)], totals,
                    args.min_cpm, args.min_count, args.min_replicates,
                )
                row["conditions"][condition][assay] = call
                condition_table.append({
                    **metadata, "condition": condition, "assay": assay,
                    "state": call["state"], "replicates": call["replicates"],
                    "passing_replicates": call["passing_replicates"],
                    "mean_count": _format(call["mean_count"]),
                    "mean_cpm": _format(call["mean_cpm"]),
                    "sample_counts": ";".join(
                        f"{sample['sample']}={_format(sample['count'])}"
                        for sample in call["samples"]
                    ),
                    "sample_cpms": ";".join(
                        f"{sample['sample']}={_format(sample['cpm']) or 'NA'}"
                        for sample in call["samples"]
                    ),
                })
        for contrast in contrasts:
            key = (feature_id, contrast)
            result = effects["deltate"].get(key)
            row["contrasts"][contrast] = {}
            xtail = effects["xtail"].get(key)
            riborex = effects["riborex"].get(key)
            for assay in EFFECT_ASSAYS:
                call = effect_call(result, assay, args.padj_cutoff, args.log2fc_cutoff)
                if assay == "TE":
                    call["xtail_te_log2fc"] = optional_float(
                        xtail.get("log2FC_TE_final") if xtail else None
                    )
                    call["xtail_te_padj"] = optional_float(
                        xtail.get("pvalue_adjusted") if xtail else None
                    )
                    call["riborex_te_log2fc"] = optional_float(
                        riborex.get("log2FC") if riborex else None
                    )
                    call["riborex_te_padj"] = optional_float(
                        riborex.get("pvalue_adjusted") if riborex else None
                    )
                row["contrasts"][contrast][assay] = call
                contrast_table.append({
                    **metadata, "contrast": contrast, "assay": assay,
                    "state": call["state"], "log2fc": _format(call["log2fc"]),
                    "padj": _format(call["padj"]), "method": "deltaTE",
                    "xtail_te_log2fc": _format(call.get("xtail_te_log2fc")),
                    "xtail_te_padj": _format(call.get("xtail_te_padj")),
                    "riborex_te_log2fc": _format(call.get("riborex_te_log2fc")),
                    "riborex_te_padj": _format(call.get("riborex_te_padj")),
                })
        rows.append(row)
    return rows, ordered_conditions, condition_table, contrast_table


def write_browser_table(path, manifest):
    tracks = []
    for track in manifest["tracks"]:
        kind = track["kind"]
        tracks.append({
            "path": f"browser/{track['file']}", "kind": kind,
            "condition": track["label"] if kind == "condition" else "",
            "contrast": track["label"] if kind == "contrast" else "",
            "assay": track["assay"],
            "direction": track["state"], "feature_count": track["feature_count"],
        })
    _write_tsv(path, TRACK_FIELDS, tracks)


def write_excel(path, rows, coordinates, conditions, contrasts, args):
    """Write the five HTML matrix views as native, sortable Excel sheets."""
    if len(rows) > 1_048_575:
        raise ValueError(
            "The cross-condition matrix exceeds Excel's 1,048,575 data-row limit; "
            "use the TSV files instead."
        )
    if max(len(conditions), len(contrasts)) > 16_382:
        raise ValueError(
            "The cross-condition matrix exceeds Excel's 16,382 data-column limit; "
            "use the TSV files instead."
        )
    workbook = xlsxwriter.Workbook(path, {"constant_memory": True})
    header = workbook.add_format({
        "bold": True, "font_color": "#FFFFFF", "bg_color": "#142744",
        "border": 1, "border_color": "#D9E2ED", "align": "center",
    })
    text_header = workbook.add_format({
        "bold": True, "font_color": "#FFFFFF", "bg_color": "#142744",
        "border": 1, "border_color": "#D9E2ED", "align": "left",
    })
    state_formats = {
        "detected": workbook.add_format({"bg_color": "#BCE8D3", "font_color": "#06422A", "align": "center"}),
        "not_detected": workbook.add_format({"bg_color": "#E5EAF0", "font_color": "#4B5563", "align": "center"}),
        "uncertain": workbook.add_format({"bg_color": "#FFF0BC", "font_color": "#6B4700", "align": "center"}),
        "up": workbook.add_format({
            "bg_color": "#F7C8C4", "font_color": "#8C1919", "align": "center",
            "num_format": '"↑ "0.0;"↓ "0.0',
        }),
        "down": workbook.add_format({
            "bg_color": "#C7DCFA", "font_color": "#123F88", "align": "center",
            "num_format": '"↑ "0.0;"↓ "0.0',
        }),
        "not_significant": workbook.add_format({"bg_color": "#F1F3F5", "font_color": "#56606A", "align": "center"}),
        "not_tested": workbook.add_format({"bg_color": "#E7EAF0", "font_color": "#56606A", "align": "center"}),
    }
    labels = {
        "detected": "Detected", "not_detected": "Not detected",
        "uncertain": "Uncertain", "not_significant": "No directional call",
        "not_tested": "Not tested",
    }

    sheet_specs = [
        ("RNA_detection", "condition", "RNA", conditions),
        ("RIBO_detection", "condition", "RIBO", conditions),
        ("RNA_change", "contrast", "RNA", contrasts),
        ("RIBO_change", "contrast", "RIBO", contrasts),
        ("TE_change", "contrast", "TE", contrasts),
    ]
    for sheet_name, kind, assay, columns in sheet_specs:
        worksheet = workbook.add_worksheet(sheet_name)
        worksheet.freeze_panes(1, 2)
        worksheet.set_column(0, 0, 20)
        worksheet.set_column(1, 1, 38)
        if columns:
            worksheet.set_column(2, len(columns) + 1, 20)
        worksheet.write_string(0, 0, "Locus tag", text_header)
        worksheet.write_string(0, 1, "Identifier", text_header)
        for column, value in enumerate(columns, 2):
            worksheet.write_string(0, column, value, header)
        for row_number, row in enumerate(rows, 1):
            metadata = _metadata(row["feature_id"], coordinates)
            # write_string prevents identity values beginning with '=' from being
            # interpreted as formulas by spreadsheet applications.
            worksheet.write_string(row_number, 0, metadata["locus_tag"])
            worksheet.write_string(row_number, 1, row["feature_id"])
            calls = row["conditions" if kind == "condition" else "contrasts"]
            for offset, label in enumerate(columns, 2):
                call = calls[label][assay]
                state = call["state"]
                if state in {"up", "down"}:
                    worksheet.write_number(
                        row_number, offset, call["log2fc"], state_formats[state]
                    )
                else:
                    worksheet.write_string(
                        row_number, offset, labels[state], state_formats[state]
                    )
        last_row = max(len(rows), 1)
        last_column = len(columns) + 1
        worksheet.autofilter(0, 0, last_row, last_column)

    readme = workbook.add_worksheet("README")
    readme.set_column(0, 0, 25)
    readme.set_column(1, 1, 105)
    readme.write_row(0, 0, ["Item", "Meaning"], text_header)
    guidance = [
        ("Workbook", "The five visual sheets mirror the selectable tables in condition_overview.html. Excel filters and native sorting are enabled."),
        ("Identity columns", "Locus tag is annotation metadata (falling back to gene_id); Identifier is the exact coordinate-style key used to join counts and differential results."),
        ("Detection", f"Detected requires at least {args.min_count} counts and {args.min_cpm} CPM in at least {args.min_replicates} usable replicates. Not detected means below these thresholds at this sequencing depth; uncertain is unresolved."),
        ("Changes", f"Up/down requires adjusted p <= {args.padj_cutoff} and absolute log2FC >= {args.log2fc_cutoff}. Positive effects are higher in the left side of a contrast."),
        ("No directional call", "The result did not meet both configured cutoffs. It does not prove no biological effect."),
        ("Not tested", "No usable statistical result was available."),
        ("Detailed values", "condition_matrix.tsv and contrast_matrix.tsv contain the exact values behind these matrix views."),
    ]
    for index, (item, meaning) in enumerate(guidance, 1):
        readme.write_string(index, 0, item)
        readme.write_string(index, 1, meaning)
    readme.autofilter(0, 0, len(guidance), 1)
    readme.freeze_panes(1, 0)

    workbook.close()


def write_html(path, rows, coordinates, conditions, contrasts, args):
    """Create a self-contained, searchable report; cap DOM rows by pagination."""
    report_rows = []
    for row in rows:
        metadata = _metadata(row["feature_id"], coordinates)
        report_rows.append({
            "id": row["feature_id"], "locus_tag": metadata["locus_tag"],
            "name": metadata["name"],
            "feature_type": metadata["feature_type"],
            "aliases": coordinates.get(row["feature_id"], {}).get("aliases", []),
            "conditions": row["conditions"], "contrasts": row["contrasts"],
        })
    payload = json.dumps({
        "conditions": conditions, "contrasts": contrasts, "features": report_rows,
        "methods": {"deltate": True, "xtail": True,
                    "riborex": getattr(args, "riborex", None) is not None},
        "thresholds": {
            "min_count": args.min_count, "min_cpm": args.min_cpm,
            "min_replicates": args.min_replicates,
            "padj_cutoff": args.padj_cutoff, "log2fc_cutoff": args.log2fc_cutoff,
        },
    }, ensure_ascii=False, separators=(",", ":")).replace("<", "\\u003c")
    template = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>HRIBO cross-condition report</title>
<style>
:root{font-family:system-ui,-apple-system,Segoe UI,sans-serif;color:#17243b;background:#f5f8fc;--locus-width:180px}
body{margin:0}header{background:#142744;color:#fff;padding:1.5rem max(1rem,calc((100vw - 1400px)/2))}
h1{font-size:1.7rem;margin:.2rem 0}h2{font-size:1.2rem}p{line-height:1.5}main{max-width:1400px;margin:auto;padding:1rem}
.card{background:#fff;border:1px solid #d9e2ed;border-radius:10px;padding:1rem;margin-bottom:1rem;box-shadow:0 2px 8px #15294c0c}
.links a{margin-right:1rem}.controls{display:flex;gap:.7rem;flex-wrap:wrap;align-items:end}.controls label{display:grid;gap:.3rem;font-weight:600;font-size:.9rem}
input,select,button{font:inherit;padding:.45rem;border:1px solid #98abc2;border-radius:5px;background:#fff}button{cursor:pointer;background:#e8f0ff}
.scroller{overflow:auto;max-height:65vh}table{border-collapse:separate;border-spacing:0;width:100%}th,td{border-bottom:1px solid #dce4ee;padding:.45rem .6rem;text-align:left;white-space:nowrap}
thead th{position:sticky;top:0;background:#eaf0f7;z-index:2}.identity{position:sticky;background:#fff;z-index:1}.locus-tag{left:0;box-sizing:border-box;width:var(--locus-width);min-width:var(--locus-width);max-width:var(--locus-width);overflow:hidden;text-overflow:ellipsis}.identifier{left:var(--locus-width);min-width:250px}
thead .identity{z-index:3;background:#eaf0f7}.identity button{border:0;background:transparent;color:#145c91;text-align:left;padding:0;font-weight:600;user-select:text}.missing{color:#77808b}
.sort-button{display:flex;width:100%;gap:.4rem;align-items:center;border:0;background:transparent;color:inherit;padding:0;font-weight:700;text-align:left}
.sort-indicator{min-width:1em;color:#145c91}.sort-button:hover,.sort-button:focus-visible{color:#145c91;text-decoration:underline}
.state{text-align:center;font-weight:650;min-width:92px;border-left:2px solid #fff}.detected{background:#bce8d3;color:#06422a}.not_detected{background:#e5eaf0;color:#4b5563}
.uncertain{background:#fff0bc;color:#6b4700}.up{background:#f7c8c4;color:#8c1919}.down{background:#c7dcfa;color:#123f88}
.not_significant{background:#f1f3f5;color:#56606a}.not_tested{background:repeating-linear-gradient(45deg,#e7eaf0,#e7eaf0 6px,#f7f8fa 6px,#f7f8fa 12px);color:#56606a}
.legend span{display:inline-block;padding:.25rem .55rem;border-radius:4px;margin:.15rem}.muted{color:#52637a}.pager{display:flex;gap:.8rem;align-items:center;margin-top:.8rem}
#detail{white-space:pre-wrap;line-height:1.5;max-height:38vh;overflow:auto}.small{font-size:.88rem}
@media (max-width:600px){:root{--locus-width:120px}.identifier{position:static;left:auto;box-sizing:border-box;width:150px;min-width:150px;max-width:150px;overflow:hidden;text-overflow:ellipsis}thead .identifier{position:sticky;top:0}}
</style></head><body>
<header><h1>Cross-condition expression and translation</h1><p>Explore feature detection and changes across conditions. Detection is not proof of biological presence or absence.</p></header>
<main><section class="card"><h2>How to read this report</h2>
<p>Only features in the differential-expression count matrix are shown. A predicted ORF not counted there is outside this report, not classified as absent or unchanged.</p>
<p>RNA and RIBO detection are separate descriptive calls. CPM uses the total counts assigned to the selected features in each sample, separately by assay. A feature is detected when both the minimum raw count and CPM are met in at least the configured number of usable replicates. A zero-assigned-count sample is unusable. “Not detected” means below these thresholds at this sequencing depth; “uncertain” includes discordant or insufficient replicates.</p>
<p>Change calls compare the left condition with the right condition in each contrast. Red is higher in the left condition; blue is lower. “No directional call” means the result did not meet both configured significance and effect-size cutoffs; it does not mean unchanged. “Not tested” includes missing or NA model results. RNA, RIBO, and translation efficiency (TE) calls use deltaTE; __SUPPLEMENTARY_METHODS__</p>
<p id="thresholds" class="small"></p><p class="links"><a href="condition_overview.xlsx">Open the same views in Excel</a><a href="condition_matrix.tsv">Condition table (TSV)</a><a href="contrast_matrix.tsv">Contrast table (TSV)</a><a href="browser_tracks.tsv">Genome-browser tracks (TSV)</a></p></section>
<section class="card"><div class="controls"><label>Search locus tag, identifier, or name<input id="search" type="search" placeholder="Locus tag, name, or coordinate"></label>
<label>View<select id="view"><option value="RNA:condition">RNA detection</option><option value="RIBO:condition">RIBO detection</option><option value="RNA:contrast">RNA change</option><option value="RIBO:contrast">RIBO change</option><option value="TE:contrast">TE change</option></select></label>
<label>Show<select id="filter"><option value="all">All features</option><option value="variable">Different states across columns</option><option value="specific">Detected in exactly one condition</option><option value="changed">Up or down in any contrast</option></select></label></div>
<p id="summary" class="muted"></p><div class="legend" id="legend"></div><div class="scroller"><table><thead><tr id="columns"></tr></thead><tbody id="matrix"></tbody></table></div>
<div class="pager"><button id="previous">Previous</button><span id="page"></span><button id="next">Next</button></div></section>
<section class="card"><h2>Selected feature</h2><div id="detail" class="muted">Select an identifier in the matrix to inspect replicate counts and each method's effect estimates.</div></section>
</main><script type="application/json" id="report-data">__DATA__</script><script>
const data=JSON.parse(document.getElementById('report-data').textContent);
const search=document.getElementById('search'),view=document.getElementById('view'),filter=document.getElementById('filter');
const matrix=document.getElementById('matrix'),columns=document.getElementById('columns'),summary=document.getElementById('summary');
const detail=document.getElementById('detail'),legend=document.getElementById('legend'),pageLabel=document.getElementById('page');
let page=0;let sortKey=null;let sortDirection='asc';const pageSize=100;const labels={detected:'Detected',not_detected:'Not detected',uncertain:'Uncertain',up:'Up',down:'Down',not_significant:'No directional call',not_tested:'Not tested'};
const stateRank={condition:{not_tested:0,not_detected:1,uncertain:2,detected:3},contrast:{not_tested:0,not_significant:1,down:2,up:3}};
document.getElementById('thresholds').textContent=`Detection: ≥${data.thresholds.min_count} counts and ≥${data.thresholds.min_cpm} CPM in ≥${data.thresholds.min_replicates} replicates. Change: adjusted p ≤${data.thresholds.padj_cutoff} and |log₂FC| ≥${data.thresholds.log2fc_cutoff}.`;
function stateOf(row,column,assay,kind){return ((kind==='condition'?row.conditions:row.contrasts)[column]||{})[assay]?.state||'not_tested'}
function chooseRows(assay,kind,cols){const query=search.value.trim().toLowerCase();return data.features.filter(row=>{
 if(query&&!`${row.locus_tag} ${row.id} ${row.name} ${row.feature_type} ${row.aliases.join(' ')}`.toLowerCase().includes(query))return false;
 const states=cols.map(col=>stateOf(row,col,assay,kind));
 if(filter.value==='variable')return new Set(states).size>1;
 if(filter.value==='specific')return kind==='condition'&&states.filter(x=>x==='detected').length===1&&states.every(x=>x==='detected'||x==='not_detected');
 if(filter.value==='changed')return kind==='contrast'&&states.some(x=>x==='up'||x==='down');
 return true;
})}
function cellOf(row,column,assay,kind){return ((kind==='condition'?row.conditions:row.contrasts)[column]||{})[assay]||{}}
function textOrder(left,right){return String(left??'').localeCompare(String(right??''),'en',{numeric:true,sensitivity:'base'})}
function numericOrder(left,right,direction){const leftMissing=!Number.isFinite(left),rightMissing=!Number.isFinite(right);if(leftMissing!==rightMissing)return leftMissing?1:-1;if(leftMissing)return 0;return direction*(left-right)}
function cellSortKey(assay,kind,column){return `cell:${assay}:${kind}:${column}`}
function sortRows(rows,assay,kind,cols){const sortedColumn=cols.find(column=>sortKey===cellSortKey(assay,kind,column));const identitySort=sortKey==='locus_tag'||sortKey==='identifier';if(!identitySort&&!sortedColumn)return rows;
 const direction=sortDirection==='asc'?1:-1;return rows.map((row,index)=>({row,index})).sort((left,right)=>{let order=0;
  if(sortKey==='locus_tag'){const leftMissing=!left.row.locus_tag,rightMissing=!right.row.locus_tag;if(leftMissing!==rightMissing)order=leftMissing?1:-1;else order=direction*(textOrder(left.row.locus_tag,right.row.locus_tag)||textOrder(left.row.id,right.row.id))}
  else if(sortKey==='identifier')order=direction*textOrder(left.row.id,right.row.id);
  else{const leftCall=cellOf(left.row,sortedColumn,assay,kind),rightCall=cellOf(right.row,sortedColumn,assay,kind);order=direction*((stateRank[kind][leftCall.state||'not_tested']??-1)-(stateRank[kind][rightCall.state||'not_tested']??-1));
   if(!order){const leftValue=kind==='condition'?leftCall.mean_cpm:Number.isFinite(leftCall.log2fc)?Math.abs(leftCall.log2fc):null,rightValue=kind==='condition'?rightCall.mean_cpm:Number.isFinite(rightCall.log2fc)?Math.abs(rightCall.log2fc):null;order=numericOrder(leftValue,rightValue,direction)}}
  return order||left.index-right.index;
 }).map(item=>item.row)}
function changeSort(key){if(sortKey===key)sortDirection=sortDirection==='asc'?'desc':'asc';else{sortKey=key;sortDirection=key==='locus_tag'||key==='identifier'?'asc':'desc'}page=0;render()}
function sortHeader(label,key,className=''){const th=document.createElement('th');th.className=className;th.scope='col';const active=sortKey===key;const currentDirection=sortDirection==='asc'?'ascending':'descending';if(active)th.setAttribute('aria-sort',currentDirection);const button=document.createElement('button');button.type='button';button.className='sort-button';const nextDirection=active&&sortDirection==='asc'?'descending':'ascending';button.title=active?`Sorted ${currentDirection}; activate for ${nextDirection}`:`Sort by ${label}`;button.setAttribute('aria-label',active?`${label}, sorted ${currentDirection}. Activate for ${nextDirection}.`:`Sort by ${label}`);button.onclick=()=>changeSort(key);const text=document.createElement('span');text.textContent=label;const indicator=document.createElement('span');indicator.className='sort-indicator';indicator.setAttribute('aria-hidden','true');indicator.textContent=active?(sortDirection==='asc'?'▲':'▼'):'↕';button.append(text,indicator);th.append(button);return th}
function showFeature(row){const lines=[`Locus tag: ${row.locus_tag||'Not available'}`,`Identifier: ${row.id}`,`Display label: ${row.name||'Not available'}`,`Feature type: ${row.feature_type||'Not available'}`,''];
 for(const condition of data.conditions){lines.push(condition);for(const assay of ['RNA','RIBO']){const call=row.conditions[condition]?.[assay];if(!call)continue;
  const samples=call.samples.map(s=>`${s.sample}: ${s.count} counts, ${s.cpm===null?'NA':s.cpm.toFixed(2)} CPM`).join('; ');
  lines.push(`  ${assay}: ${labels[call.state]} — ${call.passing_replicates}/${call.replicates} usable replicates pass; mean ${call.mean_cpm===null?'NA':call.mean_cpm.toFixed(2)} CPM`);lines.push(`    ${samples||'No samples'}`);
 }}lines.push('');for(const contrast of data.contrasts){lines.push(contrast);for(const assay of ['RNA','RIBO','TE']){const call=row.contrasts[contrast]?.[assay];if(!call)continue;
  lines.push(`  ${assay} (deltaTE): ${labels[call.state]}, log₂FC ${call.log2fc??'NA'}, padj ${call.padj??'NA'}`);
  if(assay==='TE'){let supplemental=`    xTail: log₂FC ${call.xtail_te_log2fc??'NA'}, padj ${call.xtail_te_padj??'NA'}`;__RIBOREX_DETAIL__lines.push(supplemental)}
 }}detail.textContent=lines.join('\n')}
function render(){const [assay,kind]=view.value.split(':');const cols=kind==='condition'?data.conditions:data.contrasts;const selected=sortRows(chooseRows(assay,kind,cols),assay,kind,cols);const maxPage=Math.max(0,Math.ceil(selected.length/pageSize)-1);page=Math.min(page,maxPage);
 columns.replaceChildren();columns.append(sortHeader('Locus tag','locus_tag','identity locus-tag'));columns.append(sortHeader('Identifier','identifier','identity identifier'));
 for(const col of cols)columns.append(sortHeader(col,cellSortKey(assay,kind,col)));matrix.replaceChildren();
 for(const row of selected.slice(page*pageSize,(page+1)*pageSize)){const tr=document.createElement('tr');const locus=document.createElement('td');locus.className='identity locus-tag';locus.textContent=row.locus_tag||'—';locus.title=row.locus_tag||'No locus tag available';if(!row.locus_tag)locus.classList.add('missing');tr.append(locus);const identifier=document.createElement('td');identifier.className='identity identifier';const button=document.createElement('button');button.textContent=row.id;button.title=`Open details for ${row.locus_tag||row.name||row.id}`;button.onclick=()=>showFeature(row);identifier.append(button);tr.append(identifier);
  for(const col of cols){const call=(kind==='condition'?row.conditions:row.contrasts)[col]?.[assay]||{};const state=call.state||'not_tested';const td=document.createElement('td');td.className=`state ${state}`;
   if(kind==='contrast'&&(state==='up'||state==='down'))td.textContent=`${state==='up'?'↑':'↓'} ${Math.abs(call.log2fc).toFixed(1)}`;else td.textContent=labels[state];
   td.title=kind==='condition'?`${call.passing_replicates??0}/${call.replicates??0} replicates pass; mean CPM ${call.mean_cpm===null?'NA':call.mean_cpm?.toFixed(2)??'NA'}`:`log₂FC ${call.log2fc??'NA'}; padj ${call.padj??'NA'}`;tr.append(td)}matrix.append(tr)}
 summary.textContent=`${selected.length.toLocaleString()} of ${data.features.length.toLocaleString()} features shown · ${cols.length} ${kind==='condition'?'conditions':'contrasts'} · ${assay}`;
 pageLabel.textContent=`Page ${page+1} of ${maxPage+1}`;document.getElementById('previous').disabled=page===0;document.getElementById('next').disabled=page===maxPage;
 legend.replaceChildren();for(const state of kind==='condition'?['detected','not_detected','uncertain']:['up','down','not_significant','not_tested']){const span=document.createElement('span');span.className=state;span.textContent=labels[state];legend.append(span)}
}
for(const control of [search,view,filter])control.addEventListener(control===search?'input':'change',()=>{page=0;render()});
document.getElementById('previous').onclick=()=>{page--;render()};document.getElementById('next').onclick=()=>{page++;render()};render();
</script></body></html>'''
    supplementary = (
        "xTail and RiboRex TE estimates are shown in the feature details but "
        "are not combined into a new p-value."
        if getattr(args, "riborex", None) is not None
        else "xTail TE estimates are shown in the feature details but are not "
        "combined into a new p-value."
    )
    riborex_detail = (
        "supplemental+=`; RiboRex: log₂FC "
        "${call.riborex_te_log2fc??'NA'}, padj "
        "${call.riborex_te_padj??'NA'}`;"
        if getattr(args, "riborex", None) is not None
        else ""
    )
    path.write_text(
        template.replace("__SUPPLEMENTARY_METHODS__", supplementary)
        .replace("__RIBOREX_DETAIL__", riborex_detail)
        .replace("__DATA__", payload),
        encoding="utf-8",
    )


def run(args):
    counts, labels, totals = read_counts(args.counts)
    coordinates = read_annotation(args.annotation)
    effects = {
        "xtail": read_effects(args.xtail, "xtail"),
        "riborex": (
            read_effects(args.riborex, "riborex") if args.riborex else {}
        ),
        "deltate": read_effects(args.deltate, "deltate"),
    }
    contrasts = list(dict.fromkeys(args.contrasts))
    rows, conditions, condition_table, contrast_table = build_summary(
        counts, labels, totals, coordinates, contrasts, effects, args
    )
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_dir / "condition_matrix.tsv", CONDITION_FIELDS, condition_table)
    _write_tsv(output_dir / "contrast_matrix.tsv", CONTRAST_FIELDS, contrast_table)
    write_excel(
        output_dir / "condition_overview.xlsx", rows, coordinates, conditions,
        contrasts, args,
    )
    export_tracks(rows, output_dir / "browser", conditions, contrasts, coordinates)
    manifest = json.loads((output_dir / "browser" / "tracks_manifest.json").read_text())
    write_browser_table(output_dir / "browser_tracks.tsv", manifest)
    write_html(output_dir / "condition_overview.html", rows, coordinates, conditions, contrasts, args)


def main(argv=None):
    run(parse_args(argv))


if __name__ == "__main__":
    main()
