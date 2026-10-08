"""Regression evidence from the user's real TIS-WT-1 adviser report.

The enrichment-only selector retained rare 40-nt reads while discarding the
32-nt peak with about 25 times as many boundary counts. The fixture contains
all raw threeprime start-panel rows from the HTML, with no inferred fiveprime
profiles. Its source hash and relative read-total reconstruction are recorded
in JSON; the original absolute accepted-read total remains unverified.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from lib import psite


@pytest.fixture(scope="module")
def real_tis_data():
    data = json.loads((Path(__file__).parent / "fixtures" / "tis_real_profiles.json").read_text())
    coordinates = np.asarray(data["coordinates"])
    totals = {int(length): count for length, count in data["read_totals"].items()}
    profiles = {
        int(length): np.asarray(counts, dtype=float)
        for length, counts in data["profiles"].items()
    }
    geometry = psite.READ_ENDS["threeprime"]
    scores = psite.score_read_lengths(profiles, coordinates, totals, geometry, anchor="start")
    return data, coordinates, totals, profiles, geometry, scores


def test_real_rare_length_has_higher_enrichment_but_less_boundary_support(real_tis_data):
    data, _, totals, _, _, scores = real_tis_data
    by_length = {score.read_length: score for score in scores}
    assert sum(totals.values()) == data["library_total"] == 6475868
    assert set(totals) == set(range(22, 41))

    for length in (32, 40):
        score = by_length[length]
        expected = data["expected"][str(length)]
        assert score.usable
        assert (score.offset, score.anchor, score.site) == (19, "start", "P")
        assert score.peak_height == expected["peak_height"]
        assert score.background == expected["site_aligned_background"]
        assert score.sharpness == pytest.approx(expected["site_aligned_sharpness"])
        assert score.frame_fractions[0] == pytest.approx(expected["frame_zero"])
        assert score.frame_reads == expected["frame_reads"]

    assert totals[32] == 770153
    assert totals[40] == 32415
    assert by_length[32].peak_height == 120258
    assert by_length[40].peak_height == 4896
    assert by_length[40].sharpness > by_length[32].sharpness
    assert by_length[32].frame_fractions[0] > by_length[40].frame_fractions[0]


def test_real_32_and_40_pool_preserves_strong_peak_and_increases_support(real_tis_data):
    data, coordinates, totals, profiles, geometry, scores = real_tis_data
    offsets = {score.read_length: score.offset for score in scores if score.usable}
    sharpness, frames, _ = psite._evaluate_combination(
        profiles, coordinates, offsets, [32, 40], geometry, anchor="start"
    )
    pooled = psite.pool_profiles(profiles, offsets, [32, 40], geometry)

    assert sharpness == pytest.approx(601.7019230769231)
    assert float(pooled[coordinates == 0][0]) == 125154
    assert (totals[32] + totals[40]) / data["library_total"] == pytest.approx(0.12393211226664905)
    assert frames[0] == pytest.approx(0.5430369195192126)


def test_real_default_recommendation_keeps_supported_cohort_above_quality_floor(real_tis_data):
    data, coordinates, _, profiles, geometry, scores = real_tis_data
    recommendation = psite.recommend_read_lengths(scores, profiles, coordinates, geometry, anchor="start")
    expected = data["quality_floor_recommendation"]

    assert recommendation.read_lengths == [31, 32, 33, 34, 35, 36, 38, 40]
    assert recommendation.read_lengths == expected["read_lengths"]
    assert recommendation.offsets == {31: 18, 32: 19, 33: 19, 34: 19, 35: 19, 36: 19, 38: 19, 40: 19}
    assert recommendation.selection["reference_read_length"] == 32
    assert recommendation.selection["reference_sharpness"] == pytest.approx(586.6243902439024)
    assert recommendation.selection["minimum_sharpness"] == pytest.approx(293.3121951219512)
    assert recommendation.selection["min_relative_enrichment"] == 0.5
    assert recommendation.selection["strategy"] == "coverage_with_quality_floor"
    assert recommendation.sharpness == pytest.approx(expected["sharpness"])
    assert recommendation.sharpness >= recommendation.selection["minimum_sharpness"]
    assert recommendation.covered_fraction == pytest.approx(0.40585416503239413)
    assert recommendation.covered_fraction == pytest.approx(expected["covered_fraction"])
    assert recommendation.frame_reads == expected["frame_reads"]
    assert (recommendation.anchor, recommendation.site, recommendation.read_end) == ("start", "P", "threeprime")


def test_real_stricter_quality_floor_retains_only_32_and_40(real_tis_data):
    _, coordinates, _, profiles, geometry, scores = real_tis_data
    recommendation = psite.recommend_read_lengths(
        scores, profiles, coordinates, geometry, anchor="start", min_relative_enrichment=1.0
    )

    assert recommendation.selection["reference_read_length"] == 32
    assert recommendation.read_lengths == [32, 40]
    assert recommendation.sharpness == pytest.approx(601.7019230769231)
    assert recommendation.covered_fraction == pytest.approx(0.12393211226664905)
    assert recommendation.sharpness >= recommendation.selection["minimum_sharpness"]
