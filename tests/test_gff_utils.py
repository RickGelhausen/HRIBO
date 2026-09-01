"""Tests for the shared GFF attribute parser.

The cases here are the ones that broke the previous per-script implementations:
empty values, GTF2 quoting, and round-tripping the column unchanged.
"""

import pytest

import gff_utils


def test_splits_gff3_pairs_in_order():
    pairs = gff_utils.split_attributes("ID=cds1;locus_tag=b0001;Name=gene1;")
    assert pairs == [("ID", "cds1"), ("locus_tag", "b0001"), ("Name", "gene1")]


def test_keeps_empty_values():
    """Reparation writes "ORF_type=;" and the old parser lost the pairing."""
    pairs = gff_utils.split_attributes("ID=x;ORF_type=;Length=300;")
    assert pairs == [("ID", "x"), ("ORF_type", ""), ("Length", "300")]


def test_empty_value_does_not_shift_later_keys():
    parsed = gff_utils.parse_attributes("ID=x;ORF_type=;Length=300;Prob=0.9;")
    assert parsed["length"] == "300"
    assert parsed["prob"] == "0.9"


def test_parses_gtf2_quoted_form():
    pairs = gff_utils.split_attributes('gene_id "gene1"; gene_name "abc"; locus_tag "b1";')
    assert pairs == [("gene_id", "gene1"), ("gene_name", "abc"), ("locus_tag", "b1")]


def test_value_containing_an_equals_sign_survives():
    """Splitting on every "=" at once used to break these."""
    parsed = gff_utils.parse_attributes("ID=x;Note=a=b;locus_tag=b1;")
    assert parsed["note"] == "a=b"
    assert parsed["locus_tag"] == "b1"


def test_round_trip_preserves_order_and_case():
    attributes = "ID=cds1;Parent=gene1;Name=abc;"
    assert gff_utils.format_attributes(gff_utils.split_attributes(attributes)) == attributes


def test_round_trip_preserves_empty_values():
    attributes = "ID=x;ORF_type=;Prob=0.5;"
    assert gff_utils.format_attributes(gff_utils.split_attributes(attributes)) == attributes


def test_keys_are_lowercased_values_are_not():
    parsed = gff_utils.parse_attributes("ID=Xy;Name=AbC;")
    assert parsed == {"id": "Xy", "name": "AbC"}


def test_first_occurrence_of_a_repeated_key_wins():
    parsed = gff_utils.parse_attributes("locus_tag=first;locus_tag=second;")
    assert parsed["locus_tag"] == "first"


def test_first_attribute_falls_back_in_order():
    parsed = gff_utils.parse_attributes("gene_name=abc;")
    assert gff_utils.first_attribute(parsed, "name", "gene_name") == "abc"
    assert gff_utils.first_attribute(parsed, "missing", default="-") == "-"


def test_replace_attribute_is_case_insensitive_and_keeps_order():
    pairs = gff_utils.split_attributes("ID=x;Name=old;Prob=0.5;")
    updated = gff_utils.replace_attribute(pairs, "name", "new")
    assert gff_utils.format_attributes(updated) == "ID=x;Name=new;Prob=0.5;"


def test_replace_attribute_appends_when_absent():
    pairs = gff_utils.split_attributes("ID=x;")
    updated = gff_utils.replace_attribute(pairs, "Name", "abc")
    assert gff_utils.format_attributes(updated) == "ID=x;Name=abc;"


def test_normalize_gff3_attribute_keys_preserves_reserved_names():
    pairs = gff_utils.split_attributes(
        "ID=x;Name=feature;Parent=gene1;Prob=0.5;Evidence=A-1;"
    )
    assert gff_utils.normalize_gff3_attribute_keys(pairs) == [
        ("ID", "x"),
        ("Name", "feature"),
        ("Parent", "gene1"),
        ("prob", "0.5"),
        ("evidence", "A-1"),
    ]


def test_trailing_semicolon_is_optional():
    assert gff_utils.parse_attributes("ID=x;Name=y") == gff_utils.parse_attributes("ID=x;Name=y;")


def test_whitespace_around_fields_is_ignored():
    parsed = gff_utils.parse_attributes(" ID = x ; Name = y ;")
    assert parsed == {"id": "x", "name": "y"}


def test_bare_field_without_a_value_is_kept():
    assert gff_utils.split_attributes("ID=x;flag;") == [("ID", "x"), ("flag", "")]


def test_empty_attribute_column():
    assert gff_utils.split_attributes("") == []
    assert gff_utils.parse_attributes("") == {}


def test_require_attribute_exits_with_context():
    parsed = gff_utils.parse_attributes("Name=y;")
    with pytest.raises(SystemExit) as excinfo:
        gff_utils.require_attribute(parsed, "id", context="row 4")
    assert "row 4" in str(excinfo.value)
