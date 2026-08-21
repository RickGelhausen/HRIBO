"""Tests for P-site offset estimation and the TIS caller recommendation.

The profiles here are simulated with a known offset and a known set of usable
read lengths, so the tests check that the analysis recovers ground truth rather
than merely that it runs.
"""

import numpy as np
import pytest

from lib import psite

POSITIONS_OUT = 100
POSITIONS_IN = 150
COORDINATES = np.arange(-POSITIONS_OUT, POSITIONS_IN)


def make_profile(
    offset=12,
    peak=400.0,
    background=5.0,
    periodic_amplitude=40.0,
    body_length=120,
    seed=0,
):
    """A start-codon 5'-end profile for one read length.

    A ribosome with its P-site on the start codon has its 5' end at -offset, so
    that is where the initiation peak goes. Elongation signal is laid down every
    third nucleotide downstream, in the start codon's frame.
    """
    rng = np.random.default_rng(seed)
    profile = rng.poisson(background, size=COORDINATES.size).astype(float)

    profile[COORDINATES == -offset] += peak

    if periodic_amplitude:
        psites = COORDINATES + offset
        in_body = (psites >= 15) & (psites < 15 + body_length)
        in_frame = in_body & (psites % 3 == 0)
        profile[in_frame] += periodic_amplitude

    return profile


def make_flat_profile(background=5.0, seed=1):
    """A read length with no initiation signal at all."""
    rng = np.random.default_rng(seed)
    return rng.poisson(background, size=COORDINATES.size).astype(float)


# --------------------------------------------------------------------------
# Offset estimation
# --------------------------------------------------------------------------


@pytest.mark.parametrize("planted_offset", [8, 11, 12, 13, 15])
def test_offset_is_recovered(planted_offset):
    profile = make_profile(offset=planted_offset)
    estimate = psite.estimate_offset(profile, COORDINATES)
    assert estimate.offset == planted_offset
    assert estimate.sharpness > 5
    assert estimate.is_significant


def test_offset_rejected_when_implausible():
    """A peak far from the start codon is not an initiation peak."""
    profile = make_flat_profile()
    profile[COORDINATES == -60] += 500
    estimate = psite.estimate_offset(profile, COORDINATES)
    # Either no plausible offset at all, or one that fails the significance test:
    # what must not happen is a confident offset read off background noise.
    assert not estimate.is_significant


def test_offset_none_on_empty_profile():
    estimate = psite.estimate_offset(np.zeros(COORDINATES.size), COORDINATES)
    assert estimate.offset is None
    assert estimate.peak_height == 0.0
    assert not estimate.is_significant


def test_flat_profile_has_low_sharpness():
    estimate = psite.estimate_offset(make_flat_profile(), COORDINATES)
    assert not estimate.is_significant


# --------------------------------------------------------------------------
# Frame and periodicity
# --------------------------------------------------------------------------


def test_frame_fractions_find_the_planted_frame():
    profile = make_profile(offset=12, periodic_amplitude=60.0)
    fractions = psite.frame_fractions(profile, COORDINATES, 12)
    assert np.isclose(sum(fractions), 1.0)
    assert np.argmax(fractions) == 0
    assert fractions[0] > 0.5


def test_frame_fractions_flat_without_periodicity():
    profile = make_profile(offset=12, periodic_amplitude=0.0)
    fractions = psite.frame_fractions(profile, COORDINATES, 12)
    # No frame should dominate; allow slack for Poisson noise.
    assert max(fractions) < 0.45


def test_periodicity_detected_when_present():
    strong = psite.periodicity_score(
        make_profile(offset=12, periodic_amplitude=60.0), COORDINATES, 12
    )
    absent = psite.periodicity_score(
        make_profile(offset=12, periodic_amplitude=0.0), COORDINATES, 12
    )
    assert strong > absent
    assert strong > 0.3


def test_wrong_offset_shifts_the_dominant_frame():
    """An offset off by one moves the dominant frame, which the grading uses."""
    profile = make_profile(offset=12, periodic_amplitude=60.0)
    fractions = psite.frame_fractions(profile, COORDINATES, 13)
    assert np.argmax(fractions) != 0


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------


def build_library(good_lengths=(28, 29, 30), offset=12, all_lengths=range(26, 35)):
    profiles, totals = {}, {}
    for index, length in enumerate(all_lengths):
        if length in good_lengths:
            profiles[length] = make_profile(offset=offset, seed=index)
            totals[length] = 200_000
        else:
            profiles[length] = make_flat_profile(seed=100 + index)
            totals[length] = 20_000
    return profiles, totals


def test_usable_read_lengths_are_identified():
    profiles, totals = build_library()
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    usable = {s.read_length for s in scores if s.usable}
    assert usable == {28, 29, 30}


def test_offsets_are_reported_for_usable_lengths():
    profiles, totals = build_library(offset=13)
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    for score in scores:
        if score.usable:
            assert score.offset == 13


def test_unusable_lengths_carry_a_reason():
    profiles, totals = build_library()
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    for score in scores:
        if not score.usable:
            assert score.reasons, f"read length {score.read_length} rejected without a reason"


def test_low_abundance_length_is_rejected():
    profiles, totals = build_library()
    totals[28] = 1  # sharp profile, but essentially no reads
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    rejected = next(s for s in scores if s.read_length == 28)
    assert not rejected.usable
    assert any("of the library" in r for r in rejected.reasons)


# --------------------------------------------------------------------------
# Recommendation
# --------------------------------------------------------------------------


def test_recommendation_selects_the_planted_read_lengths():
    profiles, totals = build_library(good_lengths=(28, 29, 30), offset=12)
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    recommendation = psite.recommend_read_lengths(scores, profiles, COORDINATES)

    assert recommendation.has_recommendation
    assert set(recommendation.read_lengths) <= {28, 29, 30}
    assert recommendation.read_lengths
    assert all(offset == 12 for offset in recommendation.offsets.values())


def test_recommendation_pools_and_sharpens():
    """Pooling several good lengths should beat the best single length."""
    profiles, totals = build_library(good_lengths=(28, 29, 30), offset=12)
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    recommendation = psite.recommend_read_lengths(scores, profiles, COORDINATES)

    single = max(
        (s.sharpness for s in scores if s.usable and s.read_length == recommendation.read_lengths[0]),
        default=0,
    )
    if len(recommendation.read_lengths) > 1:
        assert recommendation.sharpness >= single


def test_no_recommendation_when_nothing_is_usable():
    profiles = {length: make_flat_profile(seed=length) for length in range(26, 35)}
    totals = {length: 50_000 for length in profiles}
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    recommendation = psite.recommend_read_lengths(scores, profiles, COORDINATES)

    assert not recommendation.has_recommendation
    assert recommendation.confidence == "none"
    assert recommendation.warnings
    assert any("RNA-seq" in w for w in recommendation.warnings)


def test_confidence_is_high_with_peak_and_frame():
    profiles, totals = build_library(good_lengths=(28, 29, 30), offset=12)
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    recommendation = psite.recommend_read_lengths(scores, profiles, COORDINATES)
    assert recommendation.confidence in {"high", "medium"}


def test_weak_periodicity_is_reported_not_fatal():
    """Bacterial data often lacks periodicity; that must not block a recommendation."""
    profiles, totals = {}, {}
    for index, length in enumerate(range(26, 35)):
        if length in (29, 30):
            profiles[length] = make_profile(offset=12, seed=index, periodic_amplitude=0.0)
            totals[length] = 200_000
        else:
            profiles[length] = make_flat_profile(seed=200 + index)
            totals[length] = 20_000

    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    recommendation = psite.recommend_read_lengths(scores, profiles, COORDINATES)

    assert recommendation.has_recommendation
    assert set(recommendation.read_lengths) <= {29, 30}
    assert any("periodicity" in w.lower() or "frame" in w.lower() for w in recommendation.warnings)


def test_offset_table_is_sorted():
    profiles, totals = build_library()
    scores = psite.score_read_lengths(profiles, COORDINATES, totals)
    recommendation = psite.recommend_read_lengths(scores, profiles, COORDINATES)
    table = recommendation.offset_table()
    assert table == sorted(table)


# --------------------------------------------------------------------------
# Profile shifting
# --------------------------------------------------------------------------


def test_shift_moves_the_peak_to_the_start_codon():
    profile = make_profile(offset=12)
    shifted = psite.shift_profile(profile, 12)
    assert int(np.argmax(shifted)) == int(np.where(COORDINATES == 0)[0][0])


def test_shift_preserves_length_and_is_zero_filled():
    profile = make_profile(offset=12)
    shifted = psite.shift_profile(profile, 12)
    assert shifted.size == profile.size
    assert np.all(shifted[:12] == 0)


def test_shift_by_zero_is_identity():
    profile = make_profile(offset=12)
    assert np.array_equal(psite.shift_profile(profile, 0), profile)


def test_pooling_aligns_different_offsets():
    """Read lengths with different offsets must still stack at the start codon."""
    profiles = {28: make_profile(offset=11, seed=1), 30: make_profile(offset=14, seed=2)}
    pooled = psite.pool_profiles(profiles, {28: 11, 30: 14}, [28, 30])
    assert int(np.argmax(pooled)) == int(np.where(COORDINATES == 0)[0][0])


# --------------------------------------------------------------------------
# Detector characteristics
#
# The peak test decides whether a recommendation is made at all, so its error
# rates are worth pinning down rather than checking a single lucky seed.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("background", [2, 5, 20, 100, 1000])
def test_noise_is_rarely_mistaken_for_a_peak(background):
    false_positives = 0
    trials = 200
    for seed in range(trials):
        rng = np.random.default_rng(seed)
        profile = rng.poisson(background, size=COORDINATES.size).astype(float)
        if psite.estimate_offset(profile, COORDINATES).is_significant:
            false_positives += 1
    assert false_positives / trials < 0.02, (
        f"{false_positives}/{trials} pure-noise profiles were called significant "
        f"at background {background}"
    )


@pytest.mark.parametrize("peak_multiple", [3, 5, 10, 30])
def test_real_peaks_are_detected(peak_multiple):
    background = 20
    detected = 0
    trials = 200
    for seed in range(trials):
        rng = np.random.default_rng(seed)
        profile = rng.poisson(background, size=COORDINATES.size).astype(float)
        profile[COORDINATES == -12] += background * peak_multiple
        estimate = psite.estimate_offset(profile, COORDINATES)
        if estimate.is_significant and estimate.offset == 12:
            detected += 1
    assert detected / trials > 0.95, (
        f"only {detected}/{trials} peaks at {peak_multiple}x background were detected"
    )
