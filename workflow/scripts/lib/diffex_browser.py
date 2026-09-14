"""Export cross-condition differential-expression calls as browser-ready GFF3.

Only positive detection and significant directional calls supplied by the
report are exported. An omitted feature must never be interpreted as absent or
unchanged; the report and its thresholds remain the authority for those states.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from urllib.parse import quote


ASSAYS = ("RNA", "RIBO")
CONTRAST_ASSAYS = ("RNA", "RIBO", "TE")
EFFECT_FIELDS = ("log2fc", "padj")
DETECTION_FIELDS = ("normalized_count", "replicate_fraction")


def _safe_filename_label(label: str) -> str:
    """Make a readable, collision-resistant filename token from a user label."""
    label = str(label)
    stem = re.sub(r"[^A-Za-z0-9]+", "-", label).strip("-")[:48] or "label"
    digest = hashlib.sha256(label.encode("utf-8")).hexdigest()[:8]
    return f"{stem}-{digest}"


def _escape(value: object) -> str:
    """Percent-escape a GFF3 field or attribute value, including delimiters."""
    return quote(str(value), safe="-_.~")


def _numeric_attributes(payload: dict, names: tuple[str, ...]) -> list[tuple[str, str]]:
    attributes = []
    for name in names:
        value = payload.get(name)
        if value is None or isinstance(value, bool):
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(numeric):
            attributes.append((name, str(value)))
    return attributes


def _record(feature_id: str, coordinates: dict, label: str, assay: str,
            state: str, kind: str, payload: dict) -> str:
    try:
        seqid = str(coordinates["seqid"])
        start = int(coordinates["start"])
        end = int(coordinates["end"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Invalid coordinates for feature {feature_id!r}") from error
    if not seqid or start < 1 or end < start:
        raise ValueError(f"Invalid coordinates for feature {feature_id!r}")
    strand = coordinates.get("strand", ".")
    if strand not in ("+", "-", ".", "?"):
        raise ValueError(f"Invalid strand for feature {feature_id!r}")
    feature_type = str(coordinates.get("type") or "gene")
    name = coordinates.get("name") or feature_id
    attributes = [
        ("ID", feature_id),
        ("Name", name),
        ("gene_id", feature_id),
        ("condition" if kind == "condition" else "contrast", label),
        ("assay", assay),
        ("state", state),
    ]
    attributes.extend(_numeric_attributes(
        payload, DETECTION_FIELDS if kind == "condition" else EFFECT_FIELDS
    ))
    attr_text = ";".join(f"{key}={_escape(value)}" for key, value in attributes)
    return "\t".join((
        _escape(seqid), "HRIBO", _escape(feature_type), str(start), str(end),
        ".", strand, ".", attr_text,
    ))


def export_tracks(rows, output_dir, conditions, contrasts, feature_coords) -> list[Path]:
    """Write condition and contrast GFF3 tracks plus a JSON manifest.

    ``rows`` contains one mapping per feature with ``feature_id`` and optional
    ``conditions`` / ``contrasts`` mappings. A condition maps to ``RNA`` and
    ``RIBO`` payloads whose ``state`` may be ``detected``, ``not_detected``, or
    ``uncertain``. A contrast maps to ``RNA``, ``RIBO``, and ``TE`` payloads whose
    ``state`` may be ``up``, ``down``, ``not_significant``, or ``not_tested``.
    Optional numerical payload fields are ``normalized_count`` and
    ``replicate_fraction`` for detection, and ``log2fc`` and ``padj`` for DE.

    ``feature_coords`` maps feature IDs to ``seqid``, 1-based inclusive
    ``start``/``end``, ``strand``, and optional ``type``/``name``. Features
    without coordinates are listed in the manifest and omitted from tracks.
    Every requested track is written, including tracks with no features.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    conditions = sorted({str(condition) for condition in conditions})
    contrasts = sorted({str(contrast) for contrast in contrasts})
    rows_by_id = {}
    for row in rows:
        feature_id = str(row["feature_id"])
        if feature_id in rows_by_id:
            raise ValueError(f"Duplicate feature ID: {feature_id!r}")
        rows_by_id[feature_id] = row

    track_specs = [
        ("condition", condition, assay, "detected")
        for condition in conditions for assay in ASSAYS
    ] + [
        ("contrast", contrast, assay, state)
        for contrast in contrasts for assay in CONTRAST_ASSAYS
        for state in ("up", "down")
    ]

    outputs = []
    manifest_tracks = []
    for kind, label, assay, state in track_specs:
        filename = f"{kind}__{_safe_filename_label(label)}__{assay}_{state}.gff3"
        path = output_dir / filename
        records = []
        for feature_id, row in rows_by_id.items():
            payload = row.get("conditions" if kind == "condition" else "contrasts", {})
            payload = payload.get(label, {}).get(assay, {})
            if payload.get("state") != state or feature_id not in feature_coords:
                continue
            coordinates = feature_coords[feature_id]
            records.append((
                str(coordinates["seqid"]), int(coordinates["start"]),
                int(coordinates["end"]), feature_id,
                _record(feature_id, coordinates, label, assay, state, kind, payload),
            ))
        records.sort(key=lambda entry: entry[:4])
        path.write_text(
            "##gff-version 3\n" + "".join(record[-1] + "\n" for record in records),
            encoding="utf-8",
        )
        outputs.append(path)
        manifest_tracks.append({
            "file": filename, "kind": kind, "label": label,
            "assay": assay, "state": state, "feature_count": len(records),
        })

    missing_coords = sorted(set(rows_by_id) - set(feature_coords))
    manifest_path = output_dir / "tracks_manifest.json"
    manifest_path.write_text(json.dumps({
        "format": "GFF3", "coordinate_system": "1-based inclusive",
        "interpretation": "Only detected, up, or down calls are shown; omission is not absence.",
        "tracks": manifest_tracks,
        "unmapped_feature_ids": missing_coords,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return [*outputs, manifest_path]
