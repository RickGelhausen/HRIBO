"""Counterexamples to frame-driven ranking and padded-background enrichment."""

import numpy as np
import pytest

from lib import psite


@pytest.mark.parametrize("anchor", ["start", "stop"])
@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_stronger_peak_wins_over_stronger_body_frame(anchor, read_end):
    coordinates = np.arange(-147, 150)
    geometry = psite.read_end_for_anchor(read_end, anchor)
    offset = 16 if read_end == "threeprime" else 15 if anchor == "stop" else 12
    p_positions = coordinates + geometry.sign * offset - (3 if anchor == "stop" else 0)
    body = (p_positions >= -105) & (p_positions < -15) if anchor == "stop" else (
        (p_positions >= 15) & (p_positions < 105)
    )
    profiles = {}
    for length, peak, body_amplitude in [(32, 1000, 6), (34, 250, 27)]:
        profile = np.full(coordinates.size, 5.0)
        profile[body & (p_positions % 3 == 0)] += body_amplitude
        profile[coordinates == -geometry.sign * offset] = peak
        profiles[length] = profile
    scores = psite.score_read_lengths(profiles, coordinates, {32: 160000, 34: 86000}, geometry, anchor)
    by_length = {s.read_length: s for s in scores}
    assert all(s.usable for s in scores)
    assert by_length[32].frame_bias < by_length[34].frame_bias
    assert by_length[32].score > by_length[34].score
    recommendation = psite.recommend_read_lengths(scores, profiles, coordinates, geometry, anchor)
    assert recommendation.read_lengths == [32]
    assert recommendation.sharpness == pytest.approx(200)
    assert "read length 32 is the coverage-supported reference" in recommendation.rationale[0]


def test_zero_padding_is_not_part_of_the_measured_background():
    coordinates = np.arange(-65, 61)
    profile = np.full(coordinates.size, 5.0)
    profile[coordinates == -20] = 1000
    scores = psite.score_read_lengths({32: profile}, coordinates, {32: 10000})
    assert scores[0].usable
    assert scores[0].background_positions == 16
    recommendation = psite.recommend_read_lengths(scores, {32: profile}, coordinates)
    assert recommendation.sharpness == pytest.approx(scores[0].sharpness)
    assert recommendation.sharpness == pytest.approx(200)


@pytest.mark.parametrize("anchor,offset", [("start", 12), ("stop", 16)])
def test_one_body_read_cannot_provide_high_confidence_frame_evidence(anchor, offset):
    coordinates = np.arange(-147, 150)
    profile = np.zeros(coordinates.size)
    profile[coordinates == -offset] = 1000
    # In-frame after either a P-site or A-site calibration, inside the coding body.
    p_position = -60 if anchor == "stop" else 60
    raw_position = p_position - offset + (3 if anchor == "stop" else 0)
    profile[coordinates == raw_position] = 1
    scores = psite.score_read_lengths({32: profile}, coordinates, {32: 10000}, anchor=anchor)
    assert scores[0].usable
    assert scores[0].frame_fractions == (1, 0, 0)
    recommendation = psite.recommend_read_lengths(scores, {32: profile}, coordinates, anchor=anchor)
    assert recommendation.frame_reads == 1
    assert recommendation.confidence == "medium"
    assert any("coding-body read ends" in w for w in recommendation.warnings)
