"""Regression checks for independent boundary geometry and confidence evidence."""

import numpy as np
import pytest

from lib import psite


@pytest.mark.parametrize(
    "library_type,peak,read_end,anchor,site",
    [
        ("TIS", "start", "fiveprime", "start", "P"),
        ("RIBO", "start", "fiveprime", "start", "P"),
        ("TTS", "stop", "threeprime", "stop", "A"),
    ],
)
def test_abundance_filter_retains_boundary_footprints_for_both_read_ends(
    tmp_path, library_type, peak, read_end, anchor, site,
):
    """A selected endpoint outside CDS must not erase the true boundary signal."""
    pytest.importorskip("pysam")
    pytest.importorskip("interlap")
    import simulate_riboseq as sim
    import tis_advisor as advisor

    genome = tmp_path / "genome.fa"
    annotation = tmp_path / "annotation.gff"
    bam = tmp_path / f"{library_type}-boundary-1.bam"
    sim.write_genome(genome)
    sim.write_annotation(annotation)
    sim.write_bam(
        bam, periodic=False, seed=11, peak=peak, anchor=read_end,
        body_count=0, noise=False,
    )
    args = advisor.parse_arguments([
        "-b", str(bam), "-a", str(annotation), "-g", str(genome),
        "-o", str(tmp_path / "out"), "-r", "24-35",
        "--filtering_methods", "rpkm", "length", "--rpkm_threshold", "10",
    ])
    profiles, totals, _ = advisor.build_profiles(args)
    start_axis, stop_axis = advisor.metagene_coordinates(
        args.positions_out_ORF, args.positions_in_ORF,
    )
    coordinates = stop_axis + 3 if anchor == "stop" else start_axis
    boundary_profiles = {
        end: stop if anchor == "stop" else start
        for end, (start, stop) in profiles.items()
    }
    best, comparisons = psite.compare_read_ends(
        boundary_profiles, coordinates, totals, anchor=anchor,
    )

    assert best is not None
    assert best.read_end == read_end
    assert (best.recommendation.anchor, best.recommendation.site) == (anchor, site)
    for comparison in comparisons:
        usable = {score.read_length: score.offset for score in comparison.scores if score.usable}
        assert set(usable) == sim.GOOD_LENGTHS
        if peak == "start":
            expected = {
                length: sim.PLANTED_OFFSET if comparison.read_end == "fiveprime"
                else length - 1 - sim.PLANTED_OFFSET
                for length in sim.GOOD_LENGTHS
            }
        else:
            expected = {
                length: sim.PLANTED_TTS_THREE_PRIME_A_SITE_OFFSET
                if comparison.read_end == "threeprime"
                else length - 1 - sim.PLANTED_TTS_THREE_PRIME_A_SITE_OFFSET
                for length in sim.GOOD_LENGTHS
            }
        assert usable == expected
        assert sum(profile.sum() for profile in boundary_profiles[comparison.read_end].values()) == sum(totals.values())


@pytest.mark.parametrize("anchor", ["start", "stop"])
def test_unobserved_background_cannot_support_a_boundary_recommendation(anchor):
    """A tall peak in a truncated window is not measured enrichment evidence."""
    coordinates = np.arange(-10, 11)
    profile = np.ones(coordinates.size)
    profile[coordinates == -8] = 1000

    scores = psite.score_read_lengths(
        {30: profile}, coordinates, {30: 10000}, anchor=anchor,
    )
    assert not scores[0].usable
    assert any("background" in reason for reason in scores[0].reasons)
    recommendation = psite.recommend_read_lengths(
        scores, {30: profile}, coordinates, anchor=anchor,
    )
    assert not recommendation.has_recommendation
    assert recommendation.confidence == "none"


@pytest.mark.parametrize("offset", [-60, -30, -20, 20, 30, 60])
def test_shifts_larger_than_the_profile_leave_no_observed_counts(offset):
    """An entirely displaced profile is empty regardless of shift direction."""
    profile = np.arange(20, dtype=float)
    shifted = psite.shift_profile(profile, offset)
    assert shifted.shape == profile.shape
    assert not shifted.any()


@pytest.mark.parametrize("anchor,offset", [("start", 12), ("stop", 16)])
def test_dominant_wrong_frame_does_not_upgrade_offset_confidence(anchor, offset):
    """Frame evidence contradicting the calibrated codon cannot certify its offset."""
    coordinates = np.arange(-147, 150)
    profile = np.full(coordinates.size, 5.0)
    profile[coordinates == -offset] = 1000
    psite_positions = coordinates + offset - (3 if anchor == "stop" else 0)
    if anchor == "stop":
        body = (psite_positions >= -105) & (psite_positions < -15)
    else:
        body = (psite_positions >= 15) & (psite_positions < 105)
    profile[body & (psite_positions % 3 == 1)] += 60

    scores = psite.score_read_lengths(
        {30: profile}, coordinates, {30: 10000}, anchor=anchor,
    )
    recommendation = psite.recommend_read_lengths(
        scores, {30: profile}, coordinates, anchor=anchor,
    )
    assert recommendation.has_recommendation
    assert recommendation.frame_fractions[1] > 0.8
    assert recommendation.confidence != "high"
    assert any("dominant reading frame is 1" in warning for warning in recommendation.warnings)


def test_truncated_three_prime_window_preserves_single_length_evidence_when_pooled():
    """Clipped observations must not alter enrichment or invent periodicity."""
    coordinates = np.arange(-65, 61)
    profile = np.full(coordinates.size, 5.0)
    profile[(coordinates >= -65) & (coordinates < -35)] = 50
    profile[coordinates == 20] = 1000
    geometry = psite.READ_ENDS["threeprime"]
    scores = psite.score_read_lengths({32: profile}, coordinates, {32: 10000}, geometry)
    recommendation = psite.recommend_read_lengths(
        scores, {32: profile}, coordinates, geometry,
    )

    score = scores[0]
    assert score.usable
    assert recommendation.sharpness == pytest.approx(score.sharpness)
    assert score.sharpness == pytest.approx(200)
    assert recommendation.frame_fractions == pytest.approx(score.frame_fractions)
    assert recommendation.frame_reads == score.frame_reads
    assert recommendation.periodicity == pytest.approx(score.periodicity)
    assert score.periodicity == 0


def test_zero_background_fallback_uses_the_same_observed_bins_before_and_after_pooling():
    """A sparse peak's finite enrichment must use a consistent denominator."""
    coordinates = np.arange(-100, 150)
    profile = np.zeros(coordinates.size)
    profile[coordinates == -12] = 1000
    scores = psite.score_read_lengths({32: profile}, coordinates, {32: 10000})
    recommendation = psite.recommend_read_lengths(scores, {32: profile}, coordinates)

    assert scores[0].usable
    assert scores[0].sharpness == pytest.approx(238)
    assert recommendation.sharpness == pytest.approx(scores[0].sharpness)
