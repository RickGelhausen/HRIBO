"""Regression checks from the user's real TTS-WT-1 adviser report.

The original saturated ranking selected 34 nt first, despite the much sharper
32-nt stop peak. Its greedy pool then diluted that stronger signal. Counts for
the original selected lengths are reconstructed independently from the HTML;
their provenance and the reconstruction residual are stored with the fixture.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lib import plotting, psite


@pytest.fixture(scope="module")
def real_tts_data():
    path = Path(__file__).parent / "fixtures" / "tts_real_profiles.json"
    data = json.loads(path.read_text())
    coordinates = np.asarray(data["coordinates"])
    totals = {int(length): count for length, count in data["read_totals"].items()}
    profiles = {
        length: np.zeros(coordinates.size)
        for length in totals
    }
    # Preserve the original full-range abundance denominator. Profiles absent
    # from the reconstruction provide no candidate signal in this regression.
    profiles.update({
        int(length): np.asarray(counts, dtype=float)
        for length, counts in data["profiles"].items()
    })
    return data, coordinates, totals, profiles


def test_real_stop_scores_share_the_calibrated_site_background(real_tts_data):
    data, coordinates, totals, profiles = real_tts_data
    scores = psite.score_read_lengths(
        profiles, coordinates, totals, psite.READ_ENDS["threeprime"], anchor="stop"
    )
    by_length = {score.read_length: score for score in scores}

    for length in (32, 33, 34):
        score = by_length[length]
        expected = data["expected"][str(length)]
        assert score.usable
        assert (score.offset, score.anchor, score.site) == (16, "stop", "A")
        assert score.peak_height == expected["peak_height"]
        assert score.background == expected["site_aligned_background"]
        assert score.sharpness == pytest.approx(expected["site_aligned_sharpness"])
        assert score.abundance == pytest.approx(totals[length] / data["library_total"])

    assert by_length[32].peak_height == 201591
    assert by_length[34].peak_height == 24923
    assert by_length[32].sharpness > 4 * by_length[34].sharpness


def test_real_stronger_stop_peak_is_not_diluted_by_the_greedy_seed(real_tts_data):
    data, coordinates, totals, profiles = real_tts_data
    geometry = psite.READ_ENDS["threeprime"]
    scores = psite.score_read_lengths(profiles, coordinates, totals, geometry, anchor="stop")
    recommendation = psite.recommend_read_lengths(
        scores, profiles, coordinates, geometry, anchor="stop"
    )

    assert recommendation.read_lengths == [32]
    assert recommendation.offsets == {32: 16}
    assert recommendation.sharpness == pytest.approx(198.80769230769232)
    assert recommendation.covered_fraction == pytest.approx(1294753 / 8102155)
    assert recommendation.sharpness > data["original_recommendation"]["sharpness"]
    assert "32" in recommendation.rationale[0]
    assert (recommendation.anchor, recommendation.site) == ("stop", "A")


def test_real_same_footprints_have_equal_sharpness_from_either_end(real_tts_data):
    _, coordinates, totals, profiles = real_tts_data
    threeprime = psite.score_read_lengths(
        profiles, coordinates, totals, psite.READ_ENDS["threeprime"], anchor="stop"
    )
    fiveprime_profiles = {
        length: psite.shift_profile(profile, -(length - 1))
        for length, profile in profiles.items()
    }
    fiveprime = psite.score_read_lengths(
        fiveprime_profiles, coordinates, totals, psite.READ_ENDS["fiveprime"], anchor="stop"
    )
    threeprime_by_length = {score.read_length: score for score in threeprime}
    fiveprime_by_length = {score.read_length: score for score in fiveprime}

    for length in (32, 33, 34):
        from_threeprime = threeprime_by_length[length]
        from_fiveprime = fiveprime_by_length[length]
        assert from_threeprime.offset + from_fiveprime.offset == length - 1
        assert from_fiveprime.peak_height == from_threeprime.peak_height
        assert from_fiveprime.background == from_threeprime.background
        assert from_fiveprime.sharpness == pytest.approx(from_threeprime.sharpness)
        assert from_fiveprime.frame_fractions == pytest.approx(from_threeprime.frame_fractions)


def test_real_stop_heatmap_shows_the_scored_ratio_and_raw_counts(real_tts_data):
    _, coordinates, totals, profiles = real_tts_data
    scores = psite.score_read_lengths(
        profiles, coordinates, totals, psite.READ_ENDS["threeprime"], anchor="stop"
    )
    by_length = {score.read_length: score for score in scores}
    lengths = [32, 33, 34]
    frame = pd.DataFrame({
        "coordinates": coordinates,
        **{str(length): profiles[length] for length in lengths},
    })
    figure = plotting.plot_metagene_heatmap(
        frame, frame, lengths, "Real TTS profile regression",
        offsets={length: by_length[length].offset for length in lengths},
        read_end="threeprime", anchor="stop", site="A",
        background_references={length: by_length[length].background_reference for length in lengths},
    )
    peak_column = int(np.flatnonzero(coordinates == 16)[0])
    stop_heatmap = figure.data[1]

    for row, length in enumerate(lengths):
        assert stop_heatmap.z[row][peak_column] == pytest.approx(by_length[length].sharpness)
        assert stop_heatmap.customdata[row][peak_column] == by_length[length].peak_height
    assert "diagnostic" in figure.data[0].hovertemplate
    assert "upstream site background" in stop_heatmap.hovertemplate
