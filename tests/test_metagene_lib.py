"""Regression tests for the metagene helper library.

Both behaviours covered here were broken: the read length parser rejected the
comma-and-interval syntax that the config documents, and the window equalisation
lost one side of the profile while aliasing every filled-in read length onto a
single shared array.
"""

import numpy as np
import pytest

import lib.io as io
import lib.misc as misc


@pytest.mark.parametrize(
    "spec, expected",
    [
        ("25-34", [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]),
        ("30", [30]),
        ("22,23,27,34-35", [22, 23, 27, 34, 35]),
        ("34-35,22", [22, 34, 35]),
        # Reversed intervals are accepted and normalised.
        ("35-34", [34, 35]),
        # Duplicates collapse.
        ("30,30,30", [30]),
    ],
)
def test_parse_read_lengths(spec, expected):
    assert io.parse_read_lengths(spec) == expected


def test_parse_read_lengths_returns_integers():
    """Mixing str and int previously made the result unsortable."""
    assert all(isinstance(value, int) for value in io.parse_read_lengths("22,23,34-35"))


def test_equalize_dictionary_keys_balances_both_sides():
    start = {"a": {30: np.ones(10, dtype=np.intp)}}
    stop = {"b": {31: np.ones(10, dtype=np.intp)}}

    start, stop = misc.equalize_dictionary_keys(start, stop, 4, 6)

    assert set(start) == set(stop) == {"a", "b"}
    for coverage in (start, stop):
        for chromosome in coverage:
            assert set(coverage[chromosome]) == {30, 31}


def test_equalize_dictionary_keys_preserves_existing_data():
    start = {"a": {30: np.ones(10, dtype=np.intp)}}
    stop = {"b": {31: np.ones(10, dtype=np.intp)}}

    start, stop = misc.equalize_dictionary_keys(start, stop, 4, 6)

    assert start["a"][30][0] == 1
    assert stop["b"][31][0] == 1


def test_equalize_dictionary_keys_does_not_alias_windows():
    start = {"a": {30: np.ones(10, dtype=np.intp)}}
    stop = {"b": {31: np.ones(10, dtype=np.intp)}}

    start, stop = misc.equalize_dictionary_keys(start, stop, 4, 6)
    start["a"][31][0] = 99

    assert stop["a"][31][0] == 0, "filled windows are shared across dictionaries"
    assert start["b"][31][0] == 0, "filled windows are shared across chromosomes"
    assert start["b"][30][0] == 0, "filled windows are shared across read lengths"


def test_equalize_dictionary_keys_window_length():
    start = {"a": {30: np.ones(10, dtype=np.intp)}}
    stop = {"b": {31: np.ones(10, dtype=np.intp)}}

    start, stop = misc.equalize_dictionary_keys(start, stop, 4, 6)

    assert len(stop["a"][31]) == 4 + 6
