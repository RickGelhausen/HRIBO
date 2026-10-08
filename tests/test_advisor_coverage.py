"""Boundary read-length advice retains usable library coverage at clear peaks."""

import itertools

import numpy as np
import pytest

from lib import psite


COORDINATES = np.arange(-147, 150)


def boundary_profiles(anchor, read_end, specifications):
    """Plant independently scaled peaks over nonperiodic coding-body counts.

    Specifications map length to (library reads, peak/background ratio,
    background counts). Counts in each metagene remain below that length's
    complete-library total; abundance is independent of peak enrichment.
    """
    geometry = psite.read_end_for_anchor(read_end, anchor)
    if read_end == "fiveprime":
        offset = 15 if anchor == "stop" else 12
    else:
        offset = 13 if anchor == "stop" else 16
    profiles, totals = {}, {}
    for length, (reads, ratio, background) in specifications.items():
        profile = np.full(COORDINATES.size, background, dtype=float)
        profile[COORDINATES == -geometry.sign * offset] = ratio * background
        assert profile.sum() <= reads
        profiles[length] = profile
        totals[length] = reads
    scores = psite.score_read_lengths(profiles, COORDINATES, totals, geometry, anchor)
    return geometry, offset, profiles, totals, scores


def advise(anchor, geometry, profiles, scores):
    return psite.recommend_read_lengths(
        scores, profiles, COORDINATES, geometry, anchor
    )


def similar_lengths(anchor, read_end):
    return boundary_profiles(anchor, read_end, {
        32: (30_000, 200, 10),
        33: (20_000, 195, 10),
        34: (10_000, 190, 10),
        # An abundant, independently significant but substantially weaker peak.
        31: (40_000, 20, 40),
    })


@pytest.mark.parametrize("anchor", ["start", "stop"])
@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_rare_extreme_enrichment_does_not_displace_abundant_clear_signal(anchor, read_end):
    geometry, offset, profiles, totals, scores = boundary_profiles(anchor, read_end, {
        32: (11_900, 200, 10),
        40: (500, 1000, 0.25),
        31: (87_600, 1, 10),
    })
    by_length = {score.read_length: score for score in scores}
    assert by_length[32].usable and by_length[40].usable
    assert by_length[32].abundance == pytest.approx(0.119)
    assert by_length[40].abundance == pytest.approx(0.005)
    assert by_length[40].sharpness > by_length[32].sharpness
    assert not by_length[31].usable

    result = advise(anchor, geometry, profiles, scores)
    assert result.read_lengths == [32, 40]
    assert result.offsets == {32: offset, 40: offset}
    assert result.covered_fraction == pytest.approx(12_400 / sum(totals.values()))
    assert result.selection["reference_read_length"] == 32
    assert result.selection["reference_sharpness"] == pytest.approx(200)
    assert result.selection["minimum_sharpness"] == pytest.approx(100)
    assert result.selection["reference_abundance"] == pytest.approx(0.119)
    assert result.selection["strategy"] == "coverage_with_quality_floor"
    assert (result.anchor, result.site, result.read_end) == (
        anchor, "A" if anchor == "stop" else "P", read_end,
    )


@pytest.mark.parametrize("anchor", ["start", "stop"])
@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_similar_strong_lengths_increase_coverage_without_requiring_peak_improvement(
    anchor, read_end,
):
    geometry, offset, profiles, totals, scores = similar_lengths(anchor, read_end)
    assert all(score.usable for score in scores)
    result = advise(anchor, geometry, profiles, scores)

    assert result.read_lengths == [32, 33, 34]
    assert result.offsets == {32: offset, 33: offset, 34: offset}
    assert result.covered_fraction == pytest.approx(60_000 / sum(totals.values()))
    assert result.sharpness == pytest.approx(195)
    assert result.sharpness < result.selection["reference_sharpness"]
    assert result.selection["reference_read_length"] == 32
    assert result.selection["min_relative_enrichment"] == 0.5
    # Flat bacterial coding bodies retain their usable boundary calibration.
    assert result.periodicity == 0
    assert result.frame_fractions == pytest.approx((1 / 3, 1 / 3, 1 / 3))
    assert result.confidence == "medium"
    assert any("bacterial" in warning.lower() for warning in result.warnings)


@pytest.mark.parametrize("anchor", ["start", "stop"])
@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_abundance_does_not_rescue_a_candidate_below_the_individual_quality_floor(
    anchor, read_end,
):
    geometry, _, profiles, _, scores = similar_lengths(anchor, read_end)
    by_length = {score.read_length: score for score in scores}
    result = advise(anchor, geometry, profiles, scores)
    weak = by_length[31]

    assert weak.usable
    assert weak.total_reads > max(by_length[length].total_reads for length in (32, 33, 34))
    assert weak.sharpness < result.selection["minimum_sharpness"]
    assert 31 not in result.read_lengths
    assert result.sharpness >= result.selection["minimum_sharpness"]


@pytest.mark.parametrize("anchor", ["start", "stop"])
@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_individually_eligible_lengths_must_also_preserve_the_pooled_quality_floor(
    anchor, read_end,
):
    geometry, offset, profiles, totals, _ = boundary_profiles(anchor, read_end, {
        32: (30_000, 200, 1), 33: (10_000, 140, 1),
    })
    site_coordinates = COORDINATES + geometry.sign * offset
    background_end = -40 if anchor == "stop" else -30
    background_indices = np.flatnonzero(
        (site_coordinates >= -100) & (site_coordinates <= background_end)
    )
    half = (len(background_indices) - 1) // 2
    # Each length's measured median remains 1, but disjoint background blocks
    # raise the median of their pool to 41. Single-length ratios do not predict
    # pooled enrichment when their background shapes differ.
    profiles[32][background_indices[:half]] = 40
    profiles[33][background_indices[-half:]] = 40
    scores = psite.score_read_lengths(profiles, COORDINATES, totals, geometry, anchor)
    assert all(score.usable for score in scores)
    result = advise(anchor, geometry, profiles, scores)

    assert next(score for score in scores if score.read_length == 33).sharpness >= (
        result.selection["minimum_sharpness"]
    )
    pooled = psite.pool_profiles(profiles, {32: offset, 33: offset}, [32, 33], geometry)
    at_site = pooled[COORDINATES == 0].item()
    pooled_background = np.median(
        pooled[(COORDINATES >= -100) & (COORDINATES <= background_end)]
    )
    assert at_site / pooled_background < result.selection["minimum_sharpness"]
    assert result.read_lengths == [32]
    assert result.covered_fraction == pytest.approx(0.75)


@pytest.mark.parametrize("anchor", ["start", "stop"])
@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_explicit_quality_floor_controls_the_coverage_tradeoff(anchor, read_end):
    geometry, _, profiles, _, scores = similar_lengths(anchor, read_end)
    strict = psite.recommend_read_lengths(
        scores, profiles, COORDINATES, geometry, anchor, min_relative_enrichment=1,
    )
    permissive = psite.recommend_read_lengths(
        scores, profiles, COORDINATES, geometry, anchor, min_relative_enrichment=0,
    )

    assert strict.read_lengths == [32]
    assert strict.covered_fraction == pytest.approx(0.3)
    assert strict.selection["minimum_sharpness"] == pytest.approx(200)
    assert permissive.read_lengths == [31, 32, 33, 34]
    assert permissive.covered_fraction == 1
    assert permissive.sharpness == pytest.approx(95)
    assert permissive.selection["minimum_sharpness"] == psite.MIN_PEAK_RATIO


@pytest.mark.parametrize("anchor", ["start", "stop"])
@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_recommendation_is_invariant_to_score_and_profile_insertion_order(anchor, read_end):
    geometry, _, profiles, _, scores = similar_lengths(anchor, read_end)
    expected = advise(anchor, geometry, profiles, scores)
    by_length = {score.read_length: score for score in scores}

    for order in itertools.permutations(profiles):
        permuted_profiles = {length: profiles[length] for length in order}
        permuted_scores = [by_length[length] for length in reversed(order)]
        result = advise(anchor, geometry, permuted_profiles, permuted_scores)
        assert result.read_lengths == expected.read_lengths
        assert result.offsets == expected.offsets
        assert result.selection == expected.selection
        assert result.covered_fraction == expected.covered_fraction
        assert result.sharpness == pytest.approx(expected.sharpness)
        assert result.confidence == expected.confidence


@pytest.mark.parametrize("minimum", [-0.01, 1.01, float("nan"), float("inf")])
def test_invalid_quality_fractions_are_rejected(minimum):
    geometry, _, profiles, _, scores = similar_lengths("start", "fiveprime")
    with pytest.raises(ValueError):
        psite.recommend_read_lengths(
            scores, profiles, COORDINATES, geometry, "start",
            min_relative_enrichment=minimum,
        )


@pytest.mark.parametrize("anchor,offset", [("start", 12), ("stop", 15)])
def test_missing_calibrated_zero_cannot_produce_a_recommendation(anchor, offset):
    coordinates = np.arange(-147, 0)
    profile = np.full(coordinates.size, 10.0)
    profile[coordinates == -offset] = 2000
    geometry = psite.read_end_for_anchor("fiveprime", anchor)
    scores = psite.score_read_lengths(
        {32: profile}, coordinates, {32: 10_000}, geometry, anchor,
    )
    # The raw peak and its background are observable; its corrected codon is not.
    assert scores[0].usable
    result = psite.recommend_read_lengths(
        scores, {32: profile}, coordinates, geometry, anchor,
    )
    assert not result.has_recommendation
    assert result.read_lengths == [] and result.offsets == {}
    assert result.confidence == "none"
    assert any("position zero" in warning for warning in result.warnings)
