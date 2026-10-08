"""ORFBounder JSON contracts and offset geometry, without calling any ORFs."""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from lib.orfbounder import export_orfbounder_inputs


REPO = Path(__file__).resolve().parents[1]
EXPORT_SCRIPT = REPO / "workflow/scripts/export_orfbounder_inputs.py"
ORFBOUNDER = Path(os.environ.get("HRIBO_TEST_ORFBOUNDER", REPO.parents[1] / "ORFBounder"))


def recommendation(read_end, offsets, *, library_type="TIS", confidence="high"):
    """Keep the measured codon separate from the derived alternative site."""
    is_stop = library_type == "TTS"
    lengths = sorted(offsets)
    measured = {str(length): distance for length, distance in offsets.items()}
    delta = 3 if read_end == "fiveprime" else -3
    return {
        "read_end": read_end,
        "read_lengths": lengths,
        "offsets": measured,
        "p_site_offsets": {
            length: distance - delta if is_stop else distance
            for length, distance in measured.items()
        },
        "a_site_offsets": {
            length: distance if is_stop else distance + delta
            for length, distance in measured.items()
        },
        "anchor": "stop" if is_stop else "start",
        "site": "A" if is_stop else "P",
        "confidence": confidence if lengths else "none",
        "rationale": [],
        "warnings": [] if lengths else ["No usable calibrated footprint."],
    }


def payload(library="TIS-one_filtered", *, library_type="TIS",
            fiveprime=None, threeprime=None, chosen="fiveprime"):
    offsets = {
        "fiveprime": {28: 11, 30: 12} if fiveprime is None else fiveprime,
        "threeprime": {30: 17, 33: 19} if threeprime is None else threeprime,
    }
    ends = {
        end: {
            "recommendation": recommendation(
                end, values, library_type=library_type,
                confidence="high" if end == "fiveprime" else "medium",
            ),
            # Evidence without a recommendation must not become a guessed setup.
            "read_lengths": [{"read_length": 29, "offset": 12, "usable": False}],
        }
        for end, values in offsets.items()
    }
    return {
        "library": library,
        "library_type": library_type,
        "anchor": "stop" if library_type == "TTS" else "start",
        "site": "A" if library_type == "TTS" else "P",
        "chosen_read_end": chosen,
        "recommendation": copy.deepcopy(ends[chosen or "fiveprime"]["recommendation"]),
        "read_ends": ends,
    }


def load_pair(directory, read_end):
    end_dir = directory / read_end
    return (
        json.loads((end_dir / "read_lengths.json").read_text()),
        json.loads((end_dir / "offsets.json").read_text()),
    )


def files_snapshot(directory):
    return {
        str(path.relative_to(directory)): path.read_bytes()
        for path in directory.rglob("*") if path.is_file()
    }


@pytest.mark.parametrize("library_type,site,anchor", [
    ("TIS", "P", "start"), ("RIBO", "P", "start"), ("TTS", "A", "stop"),
])
def test_exports_each_ends_own_sparse_selection_and_signed_measured_offsets(
    tmp_path, library_type, site, anchor,
):
    source = payload(f"{library_type}-one_filtered", library_type=library_type)
    unchanged = copy.deepcopy(source)
    output = tmp_path / "orfbounder"
    manifest = export_orfbounder_inputs([source], output)
    sample = f"{library_type}-one"

    assert load_pair(output, "fiveprime") == (
        {sample: "28,30"}, {sample: {"28": -11, "30": -12}}
    )
    assert load_pair(output, "threeprime") == (
        {sample: "30,33"}, {sample: {"30": 17, "33": 19}}
    )
    assert source == unchanged  # HRIBO's main JSON keeps positive distances.
    assert manifest == json.loads((output / "manifest.json").read_text())
    assert manifest["target"] == "ORFBounder"
    library = manifest["libraries"][sample]
    assert library["library"] == source["library"]
    assert library["library_type"] == library_type
    assert library["chosen_read_end"] == "fiveprime"
    for read_end, confidence, lengths in (
        ("fiveprime", "high", [28, 30]), ("threeprime", "medium", [30, 33]),
    ):
        assert library["read_ends"][read_end]["confidence"] == confidence
        assert library["read_ends"][read_end]["read_lengths"] == lengths
        assert library["read_ends"][read_end]["site"] == site
        assert library["read_ends"][read_end]["anchor"] == anchor
        assert manifest["read_ends"][read_end]["exported"] is True
        assert manifest["read_ends"][read_end]["samples"] == [sample]


def test_aggregate_cli_exports_distinct_sample_entries_and_omits_missing_advice(tmp_path):
    sources = [
        payload("TIS-one_filtered", threeprime={}),
        payload("TTS-two_unique", library_type="TTS", fiveprime={},
                threeprime={31: 11}, chosen="threeprime"),
        payload("RNA-three_unique", library_type="RIBO", fiveprime={},
                threeprime={}, chosen=None),
        payload("RIBO-four_unique", library_type="RIBO", fiveprime={32: 13},
                threeprime={30: 17}),
    ]
    recommendations = []
    for index, source in enumerate(sources):
        path = tmp_path / f"recommendation{index}.json"
        path.write_text(json.dumps(source))
        recommendations.append(str(path))
    output = tmp_path / "exports"
    result = subprocess.run(
        [sys.executable, str(EXPORT_SCRIPT), "--recommendations", *recommendations,
         "--output-dir", str(output)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert load_pair(output, "fiveprime") == (
        {"TIS-one": "28,30", "RIBO-four": "32"},
        {"TIS-one": {"28": -11, "30": -12}, "RIBO-four": {"32": -13}},
    )
    assert load_pair(output, "threeprime") == (
        {"TTS-two": "31", "RIBO-four": "30"},
        {"TTS-two": {"31": 11}, "RIBO-four": {"30": 17}},
    )
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["libraries"]["TTS-two"]["chosen_read_end"] == "threeprime"
    for end in ("fiveprime", "threeprime"):
        assert manifest["libraries"]["RNA-three"]["read_ends"][end]["exported"] is False
    assert manifest["read_ends"]["fiveprime"]["samples"] == ["RIBO-four", "TIS-one"]
    assert manifest["read_ends"]["threeprime"]["samples"] == ["RIBO-four", "TTS-two"]


def test_unevaluated_read_end_is_not_filled_and_explicit_zero_distance_is_valid(tmp_path):
    source = payload(fiveprime={30: 0})
    del source["read_ends"]["threeprime"]
    output = tmp_path / "exports"
    manifest = export_orfbounder_inputs([source], output)
    assert load_pair(output, "fiveprime") == (
        {"TIS-one": "30"}, {"TIS-one": {"30": 0}}
    )
    assert not (output / "threeprime/read_lengths.json").exists()
    assert not (output / "threeprime/offsets.json").exists()
    status = manifest["libraries"]["TIS-one"]["read_ends"]["threeprime"]
    assert status["exported"] is False
    assert status["confidence"] is None
    assert status["read_lengths"] == []
    assert "not evaluated" in status["reason"]


def test_no_advice_manifest_removes_stale_pairs_without_guessing_from_evidence(tmp_path):
    output = tmp_path / "exports"
    export_orfbounder_inputs([payload()], output)
    unrelated = output / "notes.txt"
    unrelated.write_text("Preserve this user note.\n")
    source = payload(fiveprime={}, threeprime={}, chosen=None)
    manifest = export_orfbounder_inputs([source], output)

    assert (output / "manifest.json").is_file()
    assert not list(output.glob("*/read_lengths.json"))
    assert not list(output.glob("*/offsets.json"))
    assert unrelated.read_text() == "Preserve this user note.\n"
    for end in ("fiveprime", "threeprime"):
        assert manifest["read_ends"][end]["exported"] is False
        assert manifest["read_ends"][end]["samples"] == []
        assert manifest["read_ends"][end]["read_lengths_file"] is None
        assert manifest["read_ends"][end]["offsets_file"] is None


@pytest.mark.parametrize("invalid", [
    "missing_offset", "boolean_offset", "float_offset", "negative_offset",
    "duplicate_length", "boolean_length", "zero_length", "wrong_site", "wrong_anchor",
])
def test_invalid_recommendations_fail_before_replacing_previous_exports(tmp_path, invalid):
    output = tmp_path / "exports"
    export_orfbounder_inputs([payload()], output)
    before = files_snapshot(output)
    source = payload()
    rec = source["read_ends"]["fiveprime"]["recommendation"]
    if invalid == "missing_offset":
        del rec["offsets"]["30"]
    elif invalid.endswith("_offset"):
        rec["offsets"]["30"] = {
            "boolean_offset": True, "float_offset": 12.5, "negative_offset": -12,
        }[invalid]
    elif invalid == "duplicate_length":
        rec["read_lengths"] = [28, 30, 30]
    elif invalid == "boolean_length":
        rec["read_lengths"] = [True, 30]
    elif invalid == "zero_length":
        rec["read_lengths"] = [0, 30]
    elif invalid == "wrong_site":
        rec["site"] = "A"
    elif invalid == "wrong_anchor":
        rec["anchor"] = "stop"
    source["recommendation"] = copy.deepcopy(rec)

    with pytest.raises(ValueError):
        export_orfbounder_inputs([source], output)
    assert files_snapshot(output) == before


def test_colliding_orfbounder_sample_keys_are_rejected_before_output_mutation(tmp_path):
    output = tmp_path / "exports"
    export_orfbounder_inputs([payload()], output)
    before = files_snapshot(output)
    with pytest.raises(ValueError, match="(?i)(duplicate|collision|collid)"):
        export_orfbounder_inputs(
            [payload("TIS-one_filtered"), payload("TIS-one_unique")], output
        )
    assert files_snapshot(output) == before


@pytest.mark.parametrize("end,signed", [("fiveprime", -12), ("threeprime", 17)])
def test_exported_signs_project_to_calibrated_codon_on_both_strands(
    tmp_path, end, signed,
):
    source = payload("TIS-roundtrip_filtered", fiveprime={30: 12}, threeprime={30: 17})
    output = tmp_path / "exports"
    export_orfbounder_inputs([source], output)

    lengths, offsets = load_pair(output, end)
    assert lengths == {"TIS-roundtrip": "30"}
    assert offsets == {"TIS-roundtrip": {"30": signed}}
    # Fixed 30-nt footprints at [100, 129] (+) and [200, 229] (-) must
    # project onto the independently planted codon coordinates 112 and 217.
    plus_endpoint, minus_endpoint = (100, 229) if end == "fiveprime" else (129, 200)
    assert plus_endpoint - offsets["TIS-roundtrip"]["30"] == 112
    assert minus_endpoint + offsets["TIS-roundtrip"]["30"] == 217

    # Also cross-check the current ORFBounder readers when its separate
    # checkout is available. Core geometry assertions above always run in CI.
    if not (ORFBOUNDER / "lib/io.py").is_file():
        return

    bam = tmp_path / "TIS-roundtrip_filtered.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr", "LN": 1000}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        for index, (start, reverse) in enumerate([(100, False), (200, True)]):
            read = pysam.AlignedSegment()
            read.query_name = f"strand{index}"
            read.query_sequence = "A" * 30
            read.query_qualities = pysam.qualitystring_to_array("I" * 30)
            read.flag = 16 if reverse else 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = "30M"
            read.set_tag("NH", 1)
            handle.write(read)
    pysam.index(str(bam))
    program = """
import json
import sys
from pathlib import Path
from lib.io import parse_read_lengths, parse_offset_json
from lib.alignment_reader import PositionReader
lengths = parse_read_lengths(Path(sys.argv[1]))
offsets = parse_offset_json(Path(sys.argv[2]))
reader = PositionReader(Path(sys.argv[3]), lengths, sys.argv[4], offsets)
positions, totals = reader.output()
print(json.dumps({"lengths": lengths, "offsets": offsets,
                  "positions": {strand: values for (_, strand), values in positions.items()},
                  "totals": totals}))
"""
    # A separate interpreter keeps ORFBounder's lib package independent of HRIBO's.
    environment = {**os.environ, "PYTHONPATH": str(ORFBOUNDER), "PYTHONDONTWRITEBYTECODE": "1"}
    result = subprocess.run(
        [sys.executable, "-c", program,
         str(output / end / "read_lengths.json"), str(output / end / "offsets.json"),
         str(bam), end],
        cwd=ORFBOUNDER, env=environment, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    parsed = json.loads(result.stdout.strip().splitlines()[-1])
    assert parsed["lengths"] == {"TIS-roundtrip": ["30"]}
    assert parsed["offsets"] == {"TIS-roundtrip": {"30": signed}}
    # Both exported read ends project onto the same calibrated codon on both strands.
    assert parsed["positions"] == {"+": {"112": 1}, "-": {"217": 1}}
    assert parsed["totals"] == {"chr": 2}
