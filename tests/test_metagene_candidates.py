"""Short-CDS cohort selection and audit reporting for metagene profiles."""

import pytest

from lib import annotation


def select(tmp_path, lines, **overrides):
    annotation_path = tmp_path / "candidates.gff"
    annotation_path.write_text("\n".join(lines) + "\n")
    options = {
        "read_intervals_dict": {},
        "total_counts_dict": {},
        "genome_length_dict": {"chr": 6000},
        "filtering_methods": [],
        "mapping_method": "fiveprime",
        "rpkm_threshold": 0,
        "overlap_distance": 50,
        "positions_out_ORF": 100,
        "positions_in_ORF": 150,
        "length_cutoff": 50,
        "return_selection": True,
    }
    options.update(overrides)
    return annotation.retrieve_annotation_positions(annotation_path, **options)


def cds(start, end, attributes="ID=candidate", strand="+", contig="chr"):
    return f"{contig}\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}"


def test_short_cds_is_retained_without_the_length_filter(tmp_path):
    starts, stops, selection = select(tmp_path, [cds(301, 330)])

    assert starts["+"]["chr"] == [(300, 302)]
    assert stops["+"]["chr"] == [(327, 329)]
    assert selection.iloc[0].to_dict() == {
        "contig": "chr",
        "start": 300,
        "end": 329,
        "strand": "+",
        "feature_id": "candidate",
        "length_nt": 30,
        "status": "retained",
        "reason": "",
    }


def test_legacy_length_filter_still_requires_the_inside_window(tmp_path):
    starts, stops, selection = select(
        tmp_path, [cds(301, 330)], filtering_methods=["length"]
    )

    assert starts == stops == {"-": {}, "+": {}}
    assert selection["reason"].tolist() == ["length"]


def test_short_cohort_upper_cutoff_is_exclusive(tmp_path):
    starts, _, selection = select(
        tmp_path,
        [
            cds(301, 599, "ID=below"),
            cds(1001, 1300, "ID=equal"),
            cds(2001, 2301, "ID=above"),
        ],
        cds_max_length=300,
    )

    assert starts["+"]["chr"] == [(300, 302)]
    assert selection["length_nt"].tolist() == [299, 300, 301]
    assert selection["status"].tolist() == ["retained", "excluded", "excluded"]
    assert selection["reason"].tolist() == ["", "cohort", "cohort"]


def test_short_cohort_overlap_filter_sees_long_cds_neighbors(tmp_path):
    starts, stops, selection = select(
        tmp_path,
        [cds(301, 360, "ID=short"), cds(361, 760, "ID=long")],
        cds_max_length=300,
        filtering_methods=["overlap"],
    )

    assert starts == stops == {"-": {}, "+": {}}
    assert selection["reason"].tolist() == ["overlap", "cohort"]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_boundary_exclusions_are_reported(tmp_path, capsys, strand):
    starts, stops, selection = select(tmp_path, [cds(20, 49, strand=strand)])

    assert starts == stops == {"-": {}, "+": {}}
    assert selection["reason"].tolist() == ["boundary"]
    assert ">>Entry removal based on boundary: 1" in capsys.readouterr().out


@pytest.mark.parametrize(
    "strand,start,end", [("+", 120, 149), ("-", 5852, 5881)]
)
def test_start_only_cohort_does_not_require_the_stop_window_to_fit(
    tmp_path, strand, start, end
):
    lines = [cds(start, end, strand=strand)]
    _, _, ordinary = select(tmp_path, lines)
    starts, _, start_only = select(tmp_path, lines, required_anchors=("start",))

    assert ordinary["reason"].tolist() == ["boundary"]
    assert start_only["status"].tolist() == ["retained"]
    expected = (start - 1, start + 1) if strand == "+" else (end - 3, end - 1)
    assert starts[strand]["chr"] == [expected]


@pytest.mark.parametrize("anchors", [(), ("body",)])
def test_required_anchors_must_name_an_existing_profile(tmp_path, anchors):
    with pytest.raises(ValueError, match="required_anchors"):
        select(tmp_path, [cds(301, 330)], required_anchors=anchors)


@pytest.mark.parametrize(
    "lines", [[], ["chr\ttest\ttRNA\t301\t330\t.\t+\t.\tID=rna"]]
)
def test_zero_cds_selection_has_fixed_columns(tmp_path, lines):
    starts, stops, selection = select(tmp_path, lines)

    assert starts == stops == {"-": {}, "+": {}}
    assert selection.empty
    assert selection.columns.tolist() == annotation.SELECTION_COLUMNS


def test_reporting_deduplicates_exact_cds_and_joins_identifiers(tmp_path):
    lines = [cds(301, 330, "ID=first"), cds(301, 330, "ID=second")]
    starts, _, selection = select(tmp_path, lines, filtering_methods=["overlap"])

    assert starts["+"]["chr"] == [(300, 302)]
    assert len(selection) == 1
    assert selection["feature_id"].tolist() == ["first,second"]

    legacy = select(tmp_path, lines, return_selection=False)
    assert len(legacy) == 2
    assert legacy[0]["+"]["chr"] == [(300, 302), (300, 302)]


@pytest.mark.parametrize(
    "attributes,expected",
    [
        ("ID=identifier;locus_tag=locus", "identifier"),
        ("locus_tag=locus", "locus"),
        (".", "chr:301-330:+"),
    ],
)
def test_feature_identifier_preference_and_coordinate_fallback(
    tmp_path, attributes, expected
):
    _, _, selection = select(tmp_path, [cds(301, 330, attributes)])

    assert selection["feature_id"].tolist() == [expected]


def test_missing_contig_is_an_audited_error(tmp_path):
    starts, stops, selection = select(tmp_path, [cds(301, 330, contig="absent")])

    assert starts == stops == {"-": {}, "+": {}}
    assert selection["reason"].tolist() == ["error"]


def test_rpkm_exclusion_is_in_the_selection_table(tmp_path):
    starts, stops, selection = select(
        tmp_path,
        [cds(301, 330)],
        filtering_methods=["rpkm"],
        total_counts_dict={"chr": 10},
    )

    assert starts == stops == {"-": {}, "+": {}}
    assert selection["reason"].tolist() == ["rpkm"]
