"""DeepRibo offset advice must be evidence-based and never change the workflow."""

from types import SimpleNamespace

import pytest

from tis_advisor import deepribo_asite_advice


def comparison(read_end, rows, confidence="high", frame_fractions=(0.65, 0.2, 0.15)):
    scores = [
        SimpleNamespace(
            read_length=length,
            offset=offset,
            total_reads=reads,
            usable=usable,
        )
        for length, offset, reads, usable in rows
    ]
    return SimpleNamespace(
        read_end=read_end,
        scores=scores,
        recommendation=SimpleNamespace(
            confidence=confidence, frame_fractions=frame_fractions
        ),
    )


def test_consistent_three_prime_lengths_suggest_a_manual_deepribo_change():
    comparisons = [
        comparison("fiveprime", [(28, 12, 600, True)]),
        comparison("threeprime", [(28, 16, 550, True), (30, 16, 350, True)]),
    ]

    advice = deepribo_asite_advice("RIBO-A-1", comparisons, 1000, 12)

    assert advice["suggested_offset"] == 13
    assert advice["current_offset"] == 12
    assert advice["applied"] is False
    assert advice["per_length_offsets"] == {"28": 13, "30": 13}
    assert advice["supporting_read_lengths"] == [28, 30]
    assert advice["agreement_fraction"] == 1.0
    assert advice["read_support_fraction"] == 0.9


def test_discordant_lengths_do_not_produce_a_single_offset():
    comparisons = [
        comparison("threeprime", [(28, 16, 450, True), (30, 18, 400, True)])
    ]

    advice = deepribo_asite_advice("RIBO-A-1", comparisons, 1000, 12)

    assert advice["suggested_offset"] is None
    assert advice["per_length_offsets"] == {"28": 13, "30": 15}
    assert "disagree" in advice["reason"]


def test_usable_lengths_must_cover_enough_of_advisor_accepted_reads():
    comparisons = [
        comparison("threeprime", [(28, 16, 150, True), (30, 16, 150, True)])
    ]

    advice = deepribo_asite_advice("RIBO-A-1", comparisons, 1000, 12)

    assert advice["suggested_offset"] is None
    assert advice["read_support_fraction"] == 0.3
    assert "too few" in advice["reason"]


def test_out_of_frame_signal_prevents_a_precise_deepribo_suggestion():
    comparisons = [
        comparison(
            "threeprime", [(28, 16, 550, True), (30, 16, 350, True)],
            frame_fractions=(0.15, 0.7, 0.15),
        )
    ]

    advice = deepribo_asite_advice("RIBO-A-1", comparisons, 1000, 12)

    assert advice["suggested_offset"] is None
    assert "not frame 0" in advice["reason"]


@pytest.mark.parametrize(
    ("read_count", "expected"), [(700, None), (850, 13)]
)
def test_one_read_length_needs_overwhelming_support(read_count, expected):
    comparisons = [comparison("threeprime", [(28, 16, read_count, True)])]

    advice = deepribo_asite_advice("RIBO-A-1", comparisons, 1000, 12)

    assert advice["suggested_offset"] == expected


@pytest.mark.parametrize(
    ("library", "comparisons", "reason"),
    [
        ("TIS-A-1", [comparison("threeprime", [(28, 16, 900, True)])], "not a RIBO"),
        ("RIBO-A-1", [comparison("fiveprime", [(28, 12, 900, True)])], "not evaluated"),
        (
            "RIBO-A-1",
            [comparison("threeprime", [(28, 16, 900, True)], confidence="low")],
            "not strong enough",
        ),
    ],
)
def test_non_applicable_or_weak_evidence_has_no_suggestion(
    library, comparisons, reason
):
    advice = deepribo_asite_advice(library, comparisons, 1000, 12)

    assert advice["suggested_offset"] is None
    assert reason in advice["reason"]
