"""Export advisor results as ORFBounder input JSON, without running ORFBounder."""

from __future__ import annotations

import json
from numbers import Integral
from pathlib import Path


READ_ENDS = ("fiveprime", "threeprime")


def recommendation_inputs(library, recommendation, read_end):
    """Convert a supported calibration to the two flat ORFBounder input maps.

    Advisor distances move downstream from a 5' end and upstream from a 3' end.
    ORFBounder subtracts its offset on + and adds it on -, so its fiveprime
    offsets are negative and its threeprime offsets positive. The directly
    calibrated site is used: P-site at starts, A-site at stops.
    """
    if read_end not in READ_ENDS:
        raise ValueError(f"Unsupported ORFBounder read end: {read_end}")
    lengths = recommendation.get("read_lengths", [])
    if not isinstance(lengths, list):
        raise ValueError(f"Recommended read lengths for {library} must be a list")
    if any(isinstance(length, bool) or not isinstance(length, Integral) or length <= 0 for length in lengths):
        raise ValueError(f"Recommended read lengths for {library} must be positive integers")
    if len(set(lengths)) != len(lengths):
        raise ValueError(f"Duplicate recommended read lengths for {library}")
    if not lengths:
        return {}, {}
    if recommendation.get("read_end", read_end) != read_end:
        raise ValueError(f"Recommendation read end does not match {read_end} for {library}")
    offsets = recommendation.get("offsets", {})
    if not isinstance(offsets, dict):
        raise ValueError(f"Recommended offsets for {library} must be a mapping")
    converted = {}
    for length in sorted(lengths):
        distance = offsets.get(str(length), offsets.get(length))
        if isinstance(distance, bool) or not isinstance(distance, Integral) or distance < 0:
            raise ValueError(f"Missing or invalid nonnegative integer offset for {library}/{length}")
        converted[str(length)] = int(distance) * (-1 if read_end == "fiveprime" else 1)
    return {library: ",".join(str(length) for length in sorted(lengths))}, {library: converted}


def build_orfbounder_inputs(payloads):
    """Build end-specific input pairs and a separate availability manifest."""
    inputs = {end: {"read_lengths": {}, "offsets": {}} for end in READ_ENDS}
    libraries = {}
    for payload in payloads:
        if not isinstance(payload, dict):
            raise ValueError("Every recommendation must be a JSON object")
        library = payload.get("library")
        if not isinstance(library, str) or not library:
            raise ValueError("Every recommendation must identify its library")
        # This is the sample-key convention used by ORFBounder's BAM readers.
        sample = library.split("_", 1)[0]
        if not sample or sample == "default":
            raise ValueError(f"Invalid or reserved ORFBounder sample key: {sample!r}")
        if sample in libraries:
            raise ValueError(f"Duplicate ORFBounder sample key: {sample}")
        library_type = payload.get("library_type", library.split("-", 1)[0])
        if library_type not in {"TIS", "RIBO", "TTS"}:
            raise ValueError(f"Unsupported library type for {library}: {library_type}")
        anchor, site = ("stop", "A") if library_type == "TTS" else ("start", "P")
        comparisons = payload.get("read_ends")
        if not isinstance(comparisons, dict):
            raise ValueError(f"Missing per-read-end recommendations for {library}")
        status = {
            "library": library,
            "library_type": library_type,
            "chosen_read_end": payload.get("chosen_read_end"),
            "read_ends": {},
        }
        libraries[sample] = status
        for end in READ_ENDS:
            comparison = comparisons.get(end)
            if comparison is None:
                status["read_ends"][end] = {
                    "exported": False, "confidence": None, "anchor": anchor,
                    "site": site, "read_lengths": [], "reason": "Read end was not evaluated.",
                }
                continue
            if not isinstance(comparison, dict):
                raise ValueError(f"Invalid read-end recommendation for {library}/{end}")
            recommendation = comparison.get("recommendation")
            if not isinstance(recommendation, dict):
                raise ValueError(f"Missing recommendation for {library}/{end}")
            if recommendation.get("anchor", anchor) != anchor or recommendation.get("site", site) != site:
                raise ValueError(f"Unexpected calibrated site for {library}/{end}; expected {anchor}/{site}")
            read_lengths, offsets = recommendation_inputs(sample, recommendation, end)
            inputs[end]["read_lengths"].update(read_lengths)
            inputs[end]["offsets"].update(offsets)
            exported = bool(read_lengths)
            status["read_ends"][end] = {
                "exported": exported,
                "confidence": recommendation.get("confidence"),
                "anchor": anchor,
                "site": site,
                "read_lengths": sorted(int(length) for length in recommendation.get("read_lengths", [])),
                "reason": "" if exported else "No supported read lengths and offsets were recommended.",
            }

    manifest = {
        "format_version": 1,
        "target": "ORFBounder",
        "offset_convention": "fiveprime: negative distance; threeprime: positive distance; first nucleotide of the calibrated codon",
        "read_ends": {
            end: {
                "mapping_method": end,
                "exported": bool(inputs[end]["read_lengths"]),
                "samples": sorted(inputs[end]["read_lengths"]),
                "read_lengths_file": f"{end}/read_lengths.json" if inputs[end]["read_lengths"] else None,
                "offsets_file": f"{end}/offsets.json" if inputs[end]["offsets"] else None,
            }
            for end in READ_ENDS
        },
        "libraries": libraries,
    }
    return inputs, manifest


def export_orfbounder_inputs(payloads, output_dir):
    """Write only supported input pairs, and always write their status manifest.

    Validate every input before changing existing exports. On a no-evidence
    rerun, remove the previous generated pair rather than leaving stale advice.
    Empty/default input maps would be misleading or rejected by ORFBounder.
    """
    inputs, manifest = build_orfbounder_inputs(payloads)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for end in READ_ENDS:
        directory = output_dir / end
        if inputs[end]["read_lengths"]:
            directory.mkdir(parents=True, exist_ok=True)
            for kind in ("read_lengths", "offsets"):
                (directory / f"{kind}.json").write_text(
                    json.dumps(inputs[end][kind], indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
        else:
            for kind in ("read_lengths", "offsets"):
                (directory / f"{kind}.json").unlink(missing_ok=True)
            if directory.is_dir() and not any(directory.iterdir()):
                directory.rmdir()
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    return manifest
