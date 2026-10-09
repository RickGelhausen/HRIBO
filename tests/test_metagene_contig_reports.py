"""Per-contig reports retain every view, local navigation, and support counts."""

from html.parser import HTMLParser
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pysam
import pytest

import metagene_profiling
from lib import io


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/metagene_profiling.py"


class Report(HTMLParser):
    def __init__(self, path):
        super().__init__(convert_charrefs=True)
        self.path = path
        self.text = path.read_text()
        self.links, self.scripts, self.tables = [], [], []
        self.link = self.table = self.row = self.cell = None
        self.feed(self.text)

    def handle_starttag(self, tag, attributes):
        attributes = dict(attributes)
        if tag == "a" and "href" in attributes:
            self.link = [attributes["href"], []]
        elif tag == "script" and "src" in attributes:
            self.scripts.append(attributes["src"])
        elif tag == "table":
            self.table = []
        elif tag == "tr" and self.table is not None:
            self.row = []
        elif tag in {"th", "td"} and self.row is not None:
            self.cell = []

    def handle_data(self, data):
        if self.link is not None:
            self.link[1].append(data)
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag == "a" and self.link is not None:
            self.links.append((self.link[0], "".join(self.link[1]).strip()))
            self.link = None
        elif tag in {"th", "td"} and self.cell is not None:
            self.row.append("".join(self.cell).strip())
            self.cell = None
        elif tag == "tr" and self.row is not None:
            self.table.append(self.row)
            self.row = None
        elif tag == "table" and self.table is not None:
            self.tables.append(self.table)
            self.table = None

    def local_targets(self):
        targets = set()
        for href, _ in self.links:
            address = urlsplit(href)
            if address.scheme or address.netloc or not address.path:
                continue
            target = (self.path.parent / unquote(address.path)).resolve()
            assert target.is_file(), f"Broken relative link {href!r} in {self.path}"
            targets.add(target)
        return targets

    def leaves(self):
        return {
            target for target in self.local_targets()
            if target.suffix == ".html" and target != self.path.resolve()
        }


def example_figure(contig, marker):
    figure = go.Figure(go.Scatter(x=[0, 1], y=[marker, marker + 1], mode="lines"))
    figure.update_layout(meta={"metagene_contig": contig})
    return (f"{contig} (per read length)", "fiveprime", figure)


@pytest.mark.parametrize("start_only", [False, True], ids=["general", "sorfs"])
def test_two_contig_reports_isolate_every_view_and_both_read_ends(tmp_path, start_only):
    start = {
        "chrA": {30: np.arange(1, 7), 31: np.full(6, 7)},
        "chrB": {30: np.arange(70, 76), 31: np.full(6, 80)},
    }
    stop = {
        contig: {length: values + 10 for length, values in lengths.items()}
        for contig, lengths in start.items()
    }
    figures = []
    for method in ("fiveprime", "threeprime"):
        figures.extend(metagene_profiling.create_metagene_figures(
            start, stop, [30, 31], tmp_path, method, "raw", 2, 4, [],
            start_only=start_only,
        ))
    assert {figure.layout.meta["metagene_contig"] for _, _, figure in figures} == {"chrA", "chrB"}
    output = tmp_path / "interactive_metagene_profiling.html"
    io.create_interactive_html(
        figures, "RIBO-example", output, "online", report_html="Library overview",
        contig_reports={contig: f"Support specific to {contig}" for contig in start},
    )
    index = Report(output)
    leaves = index.leaves()
    assert len(leaves) == 2
    assert "Library overview" in index.text
    assert "Plotly.newPlot" not in index.text
    expected_plots = 4 if start_only else 6
    for contig in start:
        leaf = next(Report(path) for path in leaves if f"Support specific to {contig}" in path.read_text())
        other = "chrB" if contig == "chrA" else "chrA"
        assert f"Support specific to {other}" not in leaf.text
        assert len(re.findall(r"Plotly\.newPlot\s*\(", leaf.text)) == expected_plots
        assert len(re.findall(f'"metagene_contig"\\s*:\\s*"{contig}"', leaf.text)) == expected_plots
        assert not re.search(f'"metagene_contig"\\s*:\\s*"{other}"', leaf.text)
        assert "fiveprime mapping" in leaf.text and "threeprime mapping" in leaf.text
        assert "overlaid read lengths" in leaf.text and "per read length" in leaf.text
        navigation = leaf.local_targets()
        assert output.resolve() in navigation
        assert leaves - {leaf.path.resolve()} <= navigation


def test_local_reports_escape_unusual_names_and_avoid_case_colliding_filenames(tmp_path):
    contigs = ["Chr/A&B", "chr/a&b"]
    output = tmp_path / "interactive_metagene_profiling.html"
    io.create_interactive_html(
        [example_figure(contig, index) for index, contig in enumerate(contigs)],
        "RIBO-example", output, "local",
    )
    index = Report(output)
    leaves = index.leaves()
    assert len(leaves) == 2
    assert len({leaf.name.casefold() for leaf in leaves}) == 2
    assert all(leaf.parent == tmp_path.resolve() for leaf in leaves)
    assert {label for _, label in index.links} >= set(contigs)
    javascript = tmp_path / "plotly.min.js"
    assert javascript.is_file() and javascript.stat().st_size > 10_000
    for path in leaves:
        leaf = Report(path)
        assert output.resolve() in leaf.local_targets()
        assert leaves - {path} <= leaf.local_targets()
        assert leaf.scripts
        for source in leaf.scripts:
            address = urlsplit(source)
            assert not address.scheme and not address.netloc
            assert (path.parent / unquote(address.path)).resolve() == javascript.resolve()


def test_rerun_removes_only_tracked_contig_leaves_and_restores_single_contig_inline(tmp_path):
    output = tmp_path / "interactive_metagene_profiling.html"
    figures = [example_figure("chrA", 1), example_figure("chrB", 10)]
    io.create_interactive_html(figures, "RIBO-example", output, "online")
    leaves = Report(output).leaves()
    assert len(leaves) == 2
    assert output.with_suffix(".contigs.json").is_file()
    unrelated = tmp_path / "interactive_metagene_profiling_contig_user-kept.html"
    unrelated.write_text("User report to retain.\n")

    io.create_interactive_html(figures[:1], "RIBO-example", output, "online")
    assert all(not leaf.exists() for leaf in leaves)
    assert unrelated.read_text() == "User report to retain.\n"
    inline = Report(output)
    assert "Plotly.newPlot" in inline.text
    assert re.search(r'"metagene_contig"\s*:\s*"chrA"', inline.text)
    assert not re.search(r'"metagene_contig"\s*:\s*"chrB"', inline.text)
    assert inline.leaves() == set()


def test_two_contig_cli_sorf_support_tables_and_library_cpm_remain_separate(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(
        ">chrA reference\n" + "A" * 4000 + "\n>chrB reference\n" + "A" * 4000 + "\n"
    )
    annotation = tmp_path / "annotation.gff"
    annotation.write_text(
        "##gff-version 3\n"
        "chrA\ttest\tCDS\t501\t560\t.\t+\t0\tID=a_short\n"
        "chrA\ttest\tCDS\t1001\t1399\t.\t+\t0\tID=a_long\n"
        "chrB\ttest\tCDS\t501\t560\t.\t+\t0\tID=b_short\n"
        "chrB\ttest\tCDS\t1001\t1239\t.\t+\t0\tID=b_short2\n"
    )
    bam = tmp_path / "RIBO-example.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrA", "LN": 4000}, {"SN": "chrB", "LN": 4000}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        for index, (reference, start, length) in enumerate([
            (0, 500, 30), (0, 1000, 30), (1, 500, 30), (1, 1000, 31), (1, 2000, 28),
        ]):
            read = pysam.AlignedSegment()
            read.query_name = f"read{index}"
            read.query_sequence = "A" * length
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            read.flag = 0
            read.reference_id = reference
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = f"{length}M"
            read.set_tag("NH", 1)
            handle.write(read)
    pysam.index(str(bam))
    output = tmp_path / "metagene"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-b", str(bam), "-a", str(annotation),
         "-g", str(genome), "-o", str(output), "-r", "30-33",
         "-m", "fiveprime", "threeprime", "-n", "cpm",
         "--positions_out_ORF", "100", "--positions_in_ORF", "150",
         "--filtering_methods", "length", "--length_cutoff", "50",
         "--sorf_max_length", "300", "--output_formats", "interactive",
         "--include_plotly_js", "online"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    short_root = output / "sorfs/cpm"
    index = Report(short_root / "interactive_metagene_profiling.html")
    leaves = index.leaves()
    assert len(leaves) == 2
    assert "Plotly.newPlot" not in index.text
    for path in leaves:
        leaf = Report(path)
        tables = [
            [dict(zip(table[0], row)) for row in table[1:]]
            for table in leaf.tables if table
        ]
        counts = next(table for table in tables if table and "input_cds" in table[0])
        support = next(table for table in tables if table and "contributing_cds" in table[0])
        contig = counts[0]["contig"]
        retained = 1 if contig == "chrA" else 2
        assert {row["contig"] for row in counts + support} == {contig}
        assert {row["mapping_method"] for row in counts} == {"fiveprime", "threeprime"}
        assert {row["retained_cds"] for row in counts} == {str(retained)}
        assert {row["anchor"] for row in support} == {"start"}
        pooled = [row for row in support if row["read_length"] == "all_selected"]
        assert len(pooled) == 2
        assert {row["contributing_cds"] for row in pooled} == {str(retained)}
        assert {row["raw_count_contributions"] for row in pooled} == {str(retained)}
        leaf.local_targets()  # Navigation and audit-table downloads all resolve.
    for method, relative30, relative31 in [("fiveprime", 0, 0), ("threeprime", 29, 30)]:
        workbook = short_root / f"{method}_readcounts_start.xlsx"
        chr_a = pd.read_excel(workbook, sheet_name="chrA").set_index("coordinates")
        chr_b = pd.read_excel(workbook, sheet_name="chrB").set_index("coordinates")
        # Five library reads include an unselected 28-nt read on the other contig.
        assert chr_a.loc[relative30, "30"] == 200_000
        assert chr_b.loc[relative30, "30"] == 200_000
        assert chr_b.loc[relative31, "31"] == 200_000
        assert not (short_root / f"{method}_readcounts_stop.xlsx").exists()


def test_cli_mixed_contig_and_no_evidence_reports_scope_mapping_specific_rpkm_counts(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chrA reference\n" + "A" * 3000 + "\n")
    annotation = tmp_path / "annotation.gff"
    annotation.write_text(
        "##gff-version 3\nchrA\ttest\tCDS\t1001\t1060\t.\t+\t0\tID=short\n"
    )
    bam = tmp_path / "RIBO-example.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrA", "LN": 3000}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        read = pysam.AlignedSegment()
        read.query_name = "overlapping_footprint"
        read.query_sequence = "A" * 30
        read.query_qualities = pysam.qualitystring_to_array("I" * 30)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 1050
        read.mapping_quality = 60
        read.cigarstring = "30M"
        read.set_tag("NH", 1)
        handle.write(read)
    pysam.index(str(bam))
    output = tmp_path / "metagene"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-b", str(bam), "-a", str(annotation),
         "-g", str(genome), "-o", str(output), "-r", "30-33",
         "-m", "fiveprime", "threeprime", "-n", "cpm",
         "--positions_out_ORF", "100", "--positions_in_ORF", "150",
         "--filtering_methods", "length", "rpkm", "--rpkm_threshold", "1",
         "--length_cutoff", "50", "--sorf_max_length", "300",
         "--output_formats", "interactive", "--include_plotly_js", "online"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    index = Report(output / "sorfs/cpm/interactive_metagene_profiling.html")
    leaves = index.leaves()
    assert len(leaves) == 2
    assert "Plotly.newPlot" not in index.text
    seen = set()
    for path in leaves:
        leaf = Report(path)
        tables = [
            [dict(zip(table[0], row)) for row in table[1:]]
            for table in leaf.tables if table
        ]
        counts = next(table for table in tables if table and "input_cds" in table[0])
        support = next(table for table in tables if table and "contributing_cds" in table[0])
        assert len(counts) == 1
        row = counts[0]
        contig = row["contig"]
        seen.add(contig)
        is_placeholder = contig == "[all contigs]"
        expected_method = "threeprime" if is_placeholder else "fiveprime"
        expected_retained = "0" if is_placeholder else "1"
        assert contig in {"chrA", "[all contigs]"}
        assert row["mapping_method"] == expected_method
        assert row["retained_cds"] == expected_retained
        assert row["excluded_rpkm"] == ("1" if is_placeholder else "0")
        assert {entry["contig"] for entry in support} == {contig}
        assert {entry["mapping_method"] for entry in support} == {expected_method}
        assert {entry["retained_cds"] for entry in support} == {expected_retained}
        pooled = next(entry for entry in support if entry["read_length"] == "all_selected")
        assert pooled["contributing_cds"] == expected_retained
        assert pooled["raw_count_contributions"] == expected_retained
        if is_placeholder:
            assert {entry["contributing_cds"] for entry in support} == {"0"}
            assert {entry["raw_count_contributions"] for entry in support} == {"0"}
            assert re.search(r'"metagene_contig"\s*:\s*"no_evidence"', leaf.text)
        leaf.local_targets()
    assert seen == {"chrA", "[all contigs]"}
