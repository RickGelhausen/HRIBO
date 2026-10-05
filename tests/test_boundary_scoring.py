"""Ground-truth tests for termination/A-site calibration and coding-body evidence."""

import numpy as np
import pytest

from lib import psite


COORDINATES = np.arange(-180, 120)


def stop_profile(read_end="fiveprime", offset=16, periodic=True, background=5.0):
    """A stop peak with independent upstream ORF and downstream decoy signals."""
    geometry = psite.READ_ENDS[read_end]
    profile = np.full(COORDINATES.size, background)
    profile[COORDINATES == -geometry.sign * offset] += 1000
    psites = COORDINATES + geometry.sign * offset - 3
    if periodic:
        body = (psites >= -105) & (psites < -15)
        profile[body & (psites % 3 == 0)] += 60
    # Large out-of-frame signal after the stop must never drive TTS frame bias.
    downstream = (psites >= 15) & (psites < 105)
    profile[downstream & (psites % 3 == 1)] += 200
    return profile


@pytest.mark.parametrize("read_end,offset", [("fiveprime", 16), ("threeprime", 13)])
def test_stop_peak_calibrates_the_a_site_from_both_read_ends(read_end, offset):
    profile = stop_profile(read_end, offset)
    geometry = psite.READ_ENDS[read_end]
    estimate = psite.estimate_offset(profile, COORDINATES, geometry, anchor="stop")
    assert estimate.offset == offset
    assert estimate.is_significant
    scores = psite.score_read_lengths({30: profile}, COORDINATES, {30: 10000}, geometry, "stop")
    assert scores[0].usable
    assert (scores[0].anchor, scores[0].site) == ("stop", "A")
    assert scores[0].offset == offset


@pytest.mark.parametrize("read_end,offset,p_offset", [("fiveprime", 16, 13), ("threeprime", 13, 16)])
def test_a_to_p_conversion_preserves_transcript_direction(read_end, offset, p_offset):
    geometry = psite.READ_ENDS[read_end]
    assert psite.psite_offset(offset, geometry, "stop") == p_offset
    assert geometry.site_shift(p_offset) == geometry.site_shift(offset) - 3
    assert psite.psite_offset(offset, geometry) == offset


@pytest.mark.parametrize("read_end,offset", [("fiveprime", 16), ("threeprime", 13)])
def test_termination_frame_and_periodicity_use_the_upstream_coding_body(read_end, offset):
    profile = stop_profile(read_end, offset)
    geometry = psite.READ_ENDS[read_end]
    fractions = psite.frame_fractions(profile, COORDINATES, offset, geometry, anchor="stop")
    assert np.argmax(fractions) == 0
    assert fractions[0] > 0.8
    assert psite.periodicity_score(profile, COORDINATES, offset, geometry, anchor="stop") > 0.9
    # The post-stop decoy dominates a start-style analysis of this profile.
    downstream = psite.frame_fractions(profile, COORDINATES, offset, geometry)
    assert np.argmax(downstream) == 1


@pytest.mark.parametrize("read_end,offset", [("fiveprime", 16), ("threeprime", 13)])
def test_stop_body_boundaries_are_applied_after_the_p_site_conversion(read_end, offset):
    geometry = psite.READ_ENDS[read_end]
    psites = COORDINATES + geometry.sign * offset - 3
    profile = np.zeros(COORDINATES.size)
    profile[psites == -106] = 1000  # before the body; excluded
    profile[psites == -17] = 50  # body frame 1; included
    profile[psites == -15] = 1000  # close to the terminal peak; excluded
    assert psite.frame_fractions(profile, COORDINATES, offset, geometry, anchor="stop") == (0, 1, 0)


def test_terminal_peak_alone_does_not_create_frame_or_periodicity_evidence():
    profile = np.full(COORDINATES.size, 5.0)
    profile[COORDINATES == -16] = 10000
    scores = psite.score_read_lengths({30: profile}, COORDINATES, anchor="stop")
    assert scores[0].usable
    assert scores[0].frame_fractions == pytest.approx((1 / 3, 1 / 3, 1 / 3))
    assert scores[0].periodicity == 0
    recommendation = psite.recommend_read_lengths(scores, {30: profile}, COORDINATES, anchor="stop")
    assert recommendation.has_recommendation
    assert recommendation.confidence == "medium"
    assert any("frame bias" in warning for warning in recommendation.warnings)


def test_stop_recommendation_reports_termination_and_a_site_metadata():
    profiles = {30: stop_profile()}
    best, comparisons = psite.compare_read_ends({"fiveprime": profiles}, COORDINATES, anchor="stop")
    assert best is not None
    recommendation = best.recommendation
    assert (recommendation.anchor, recommendation.site) == ("stop", "A")
    assert recommendation.offsets == {30: 16}
    assert recommendation.frame_fractions[0] > 0.8
    assert any("termination peak" in reason for reason in recommendation.rationale)
    assert any("A-site" in warning and "readthrough" in warning for warning in recommendation.warnings)
    assert psite.describe_end_choice(best, comparisons, "stop")
    assert not any("initiation" in warning for warning in recommendation.warnings)


def test_no_stop_signal_names_stop_annotation_and_termination():
    profiles = {30: np.full(COORDINATES.size, 5.0)}
    best, comparisons = psite.compare_read_ends(
        {"fiveprime": profiles, "threeprime": profiles}, COORDINATES, anchor="stop"
    )
    assert best is None
    assert psite.describe_end_choice(best, comparisons, "stop") == [
        "Neither read end produced a usable termination signal."
    ]
    for comparison in comparisons:
        recommendation = comparison.recommendation
        assert (recommendation.anchor, recommendation.site) == ("stop", "A")
        assert any("stop codons" in warning for warning in recommendation.warnings)
        assert not any("initiation" in warning or "start codons" in warning for warning in recommendation.warnings)


@pytest.mark.parametrize("anchor,read_length,offset", [("start", 12, 12), ("stop", 16, 16), ("stop", 15, 16)])
def test_offset_must_fall_inside_its_read(anchor, read_length, offset):
    profile = np.full(COORDINATES.size, 5.0)
    profile[COORDINATES == -offset] = 1000
    scores = psite.score_read_lengths({read_length: profile}, COORDINATES, anchor=anchor)
    assert not scores[0].usable
    assert scores[0].offset is None
    assert any("outside" in reason for reason in scores[0].reasons)
    recommendation = psite.recommend_read_lengths(scores, {read_length: profile}, COORDINATES, anchor=anchor)
    assert not recommendation.has_recommendation


@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_configured_disome_length_calibrates_the_same_terminating_ribosome(read_end):
    length = 58
    offset = length - 1 - 13 if read_end == "fiveprime" else 13
    profile = stop_profile(read_end, offset, periodic=False)
    geometry = psite.read_end_for_anchor(read_end, "stop", length)
    assert geometry.is_plausible(offset)
    scores = psite.score_read_lengths({length: profile}, COORDINATES, read_end=geometry, anchor="stop")
    assert scores[0].usable
    assert scores[0].offset == offset
    # Excluding the disome search window keeps the sharp peak out of background.
    assert scores[0].background == 5
    assert scores[0].z_score > 100


@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_start_defaults_retain_legacy_read_end_geometry(read_end):
    assert psite.read_end_for_anchor(read_end) == psite.READ_ENDS[read_end]
    assert psite.read_end_for_anchor(read_end, read_length=58) == psite.READ_ENDS[read_end]


def test_invalid_boundary_is_rejected_even_without_profiles():
    with pytest.raises(ValueError, match="anchor"):
        psite.compare_read_ends({}, COORDINATES, anchor="elongation")
    with pytest.raises(ValueError, match="anchor"):
        psite.score_read_lengths({}, COORDINATES, anchor="elongation")
