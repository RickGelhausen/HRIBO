"""Exercise separate short-CDS start profiles with actual indexed alignments."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import pysam
import pytest
import metagene_profiling


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/metagene_profiling.py"
READ_LENGTHS = ["30", "31", "32", "33"]
ALL_CONTIGS = "[all contigs]"


def test_start_only_report_plots_every_selected_read_length(tmp_path):
    lengths = list(range(22, 41))
    start = {"chr": {length: np.ones(6) for length in lengths}}
    figures = metagene_profiling.create_metagene_figures(
        start, {}, lengths, tmp_path, "fiveprime", "cpm", 2, 4, [],
        start_only=True,
    )
    assert len(figures) == 1
    figure = figures[0][2]
    assert [trace.name for trace in figure.data] == [f"{length} nt" for length in lengths]
    assert not (tmp_path / "fiveprime_readcounts_stop.xlsx").exists()


def write_inputs(tmp_path, cds_records, reads):
    """Use 1-based CDS coordinates and 0-based alignment starts."""
    genome = tmp_path / "genome.fa"
    genome.write_text(
        ">chr reference\n" + "A" * 6000 + "\n>other reference\n" + "A" * 1000 + "\n"
    )
    annotation = tmp_path / "annotation.gff"
    lines = ["##gff-version 3"]
    for identifier, start, end, strand in cds_records:
        lines.append(
            f"chr\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={identifier}"
        )
    annotation.write_text("\n".join(lines) + "\n")

    bam = tmp_path / "RIBO-test-1.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr", "LN": 6000}, {"SN": "other", "LN": 1000}],
    }
    references = {"chr": 0, "other": 1}
    with pysam.AlignmentFile(bam, "wb", header=header) as output:
        for index, (contig, start, length, strand) in enumerate(
            sorted(reads, key=lambda item: (references[item[0]], item[1]))
        ):
            read = pysam.AlignedSegment()
            read.query_name = f"read{index}"
            read.query_sequence = "A" * length
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            read.flag = 16 if strand == "-" else 0
            read.reference_id = references[contig]
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = f"{length}M"
            read.set_tag("NH", 1)
            output.write(read)
    pysam.index(str(bam))
    return bam, annotation, genome


def run_metagene(tmp_path, cds_records, reads):
    bam, annotation, genome = write_inputs(tmp_path, cds_records, reads)
    output = tmp_path / "metagene"
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "-b", str(bam), "-a", str(annotation), "-g", str(genome),
            "-o", str(output), "-r", "30-33",
            "-m", "fiveprime", "threeprime", "-n", "cpm",
            "--positions_out_ORF", "100", "--positions_in_ORF", "150",
            "--filtering_methods", "length", "--length_cutoff", "50",
            "--sorf_max_length", "300", "--output_formats", "interactive",
            "--include_plotly_js", "integrated",
        ],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return output


def report_tables(root):
    return {
        name: pd.read_csv(root / f"{name}.tsv", sep="\t", keep_default_na=False)
        for name in ("candidates", "candidate_counts", "candidate_support")
    }


def support_row(support, method, read_length):
    matches = support.loc[
        (support["contig"] == ALL_CONTIGS)
        & (support["mapping_method"] == method)
        & (support["anchor"] == "start")
        & (support["read_length"] == str(read_length))
    ]
    assert len(matches) == 1
    return matches.iloc[0]


def test_short_cds_profiles_selection_support_orientation_and_library_cpm(tmp_path):
    cds_records = [
        ("plus_short", 501, 560, "+"),
        ("minus_short", 1001, 1060, "-"),
        ("edge_short", 1601, 1660, "+"),
        ("below_cutoff", 2101, 2399, "+"),
        ("at_cutoff", 2901, 3200, "+"),
        ("long", 3801, 4199, "+"),
    ]
    reads = [
        ("chr", 500, 30, "+"),
        ("chr", 508, 31, "+"),
        ("chr", 1030, 30, "-"),
        ("chr", 1032, 33, "-"),
        # This alignment overlaps [-100,150), but only its 3' end is inside.
        ("chr", 1485, 30, "+"),
        ("chr", 2100, 32, "+"),
        ("chr", 2900, 31, "+"),
        ("chr", 3800, 30, "+"),
        # Unselected lengths still contribute to the complete-library denominator.
        ("chr", 500, 28, "+"),
        ("other", 300, 25, "+"),
    ]
    output = run_metagene(tmp_path, cds_records, reads)
    general = report_tables(output)
    short = report_tables(output / "sorfs")

    for method in ("fiveprime", "threeprime"):
        general_candidates = general["candidates"].loc[
            general["candidates"]["mapping_method"] == method
        ].set_index("feature_id")
        short_candidates = short["candidates"].loc[
            short["candidates"]["mapping_method"] == method
        ].set_index("feature_id")
        for identifier in ("plus_short", "minus_short", "edge_short"):
            assert general_candidates.loc[identifier, "reason"] == "length"
            assert short_candidates.loc[identifier, "status"] == "retained"
        assert short_candidates.loc["below_cutoff", "length_nt"] == 299
        assert short_candidates.loc["below_cutoff", "status"] == "retained"
        for identifier in ("at_cutoff", "long"):
            assert short_candidates.loc[identifier, "reason"] == "cohort"
            assert general_candidates.loc[identifier, "status"] == "retained"

        counts = short["candidate_counts"].loc[
            (short["candidate_counts"]["mapping_method"] == method)
            & (short["candidate_counts"]["contig"] == ALL_CONTIGS)
        ].iloc[0]
        assert counts["input_cds"] == 6
        assert counts["cohort_cds"] == counts["retained_cds"] == 4
        assert counts["excluded_cohort"] == 2
        assert counts["excluded_length"] == 0

        # Two lengths from each of two CDSs must not inflate the pooled CDS count.
        pooled = support_row(short["candidate_support"], method, "all_selected")
        assert pooled["retained_cds"] == 4
        assert pooled["contributing_cds"] == (3 if method == "fiveprime" else 4)
        assert pooled["raw_count_contributions"] == (5 if method == "fiveprime" else 6)
        length30 = support_row(short["candidate_support"], method, 30)
        assert length30["contributing_cds"] == (2 if method == "fiveprime" else 3)
        assert length30["raw_count_contributions"] == (2 if method == "fiveprime" else 3)
        for length in (31, 32, 33):
            support = support_row(short["candidate_support"], method, length)
            assert support["contributing_cds"] == support["raw_count_contributions"] == 1

        start = pd.read_excel(
            output / "sorfs/cpm" / f"{method}_readcounts_start.xlsx", sheet_name="chr"
        ).set_index("coordinates")
        assert start.index.tolist() == list(range(-100, 150))
        assert start.columns.tolist() == [*READ_LENGTHS, "sum"]
        assert start[READ_LENGTHS].sum().sum() == pytest.approx(
            (5 if method == "fiveprime" else 6) * 100_000
        )
        # Plus and minus 30-nt alignments both have transcript-relative ends 0/29.
        assert start.loc[0 if method == "fiveprime" else 29, "30"] == 200_000
        assert start.loc[-5 if method == "fiveprime" else 27, "33"] == 100_000
        assert start.loc[8 if method == "fiveprime" else 38, "31"] == 100_000
        assert start.loc[0 if method == "fiveprime" else 31, "32"] == 100_000
        assert start.loc[-86, "30"] == (0 if method == "fiveprime" else 100_000)
        assert not (output / "sorfs/cpm" / f"{method}_readcounts_stop.xlsx").exists()
        assert (output / "cpm" / f"{method}_readcounts_stop.xlsx").is_file()

    html = (output / "sorfs/cpm/interactive_metagene_profiling.html").read_text()
    assert "4 retained CDSs" in html
    assert "3 contributing CDSs" in html
    assert "4 contributing CDSs" in html
    assert "candidate_counts.tsv" in html and "candidate_support.tsv" in html


@pytest.mark.parametrize("cohort", ["short_without_evidence", "no_short_cds", "no_cds"])
def test_empty_short_profiles_export_explicit_support_and_html(tmp_path, cohort):
    cds_records = {
        "short_without_evidence": [("short", 501, 560, "+")],
        "no_short_cds": [("long", 501, 800, "+")],
        "no_cds": [],
    }[cohort]
    # CPM is defined, while neither a selected length nor an annotated contig has reads.
    output = run_metagene(tmp_path, cds_records, [("other", 300, 28, "+")])
    short = report_tables(output / "sorfs")
    expected_retained = int(cohort == "short_without_evidence")

    for method in ("fiveprime", "threeprime"):
        counts = short["candidate_counts"].loc[
            (short["candidate_counts"]["mapping_method"] == method)
            & (short["candidate_counts"]["contig"] == ALL_CONTIGS)
        ].iloc[0]
        assert counts["input_cds"] == len(cds_records)
        assert counts["retained_cds"] == expected_retained
        assert counts["excluded_cohort"] == int(cohort == "no_short_cds")
        for length in [*READ_LENGTHS, "all_selected"]:
            support = support_row(short["candidate_support"], method, length)
            assert support["retained_cds"] == expected_retained
            assert support["contributing_cds"] == support["raw_count_contributions"] == 0
        frame = pd.read_excel(
            output / "sorfs/cpm" / f"{method}_readcounts_start.xlsx",
            sheet_name="no_evidence",
        )
        assert frame["coordinates"].tolist() == list(range(-100, 150))
        assert frame.columns.tolist() == ["coordinates", *READ_LENGTHS, "sum"]
        assert not frame[[*READ_LENGTHS, "sum"]].to_numpy().any()

    html = (output / "sorfs/cpm/interactive_metagene_profiling.html").read_text()
    assert f"{expected_retained} retained CDSs" in html
    assert "0 contributing CDSs" in html
    assert "Candidate counts" in html and "Profile support" in html
    assert "no_evidence" in html
