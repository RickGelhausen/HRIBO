"""Regression tests for metagene settings and workflow command construction."""

import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import yaml

from lib import annotation
from lib import io as hribo_io
from lib import theme
import metagene_profiling
import tis_advisor


REPO = Path(__file__).resolve().parent.parent
RULES = (REPO / "workflow" / "rules" / "metageneprofiling.smk").read_text()
OUTPUT_FORMATS = ("interactive", "svg", "pdf", "png", "jpg")


class RecordingFigure:
    """Small Plotly stand-in that records static and HTML write requests."""

    def __init__(self):
        self.image_writes = []

    def write_image(self, path, width, height):
        self.image_writes.append((Path(path), width, height))

    def to_html(self, **kwargs):
        return '<div id="recording-figure">plot</div>'


def test_output_writer_and_schema_support_the_same_formats():
    schema = yaml.safe_load(
        (REPO / "workflow" / "schemas" / "config.schema.yaml").read_text()
    )
    advertised = schema["properties"]["metageneSettings"]["properties"][
        "outputFormats"
    ]["items"]["enum"]

    assert set(advertised) == set(OUTPUT_FORMATS)


@pytest.mark.parametrize("output_format", OUTPUT_FORMATS)
def test_every_advertised_metagene_output_format_is_written(output_format, tmp_path):
    figure = RecordingFigure()
    figures = [("chr1", "fiveprime", figure)]

    hribo_io.write_plots_to_file(
        figures,
        [output_format],
        "online",
        "RIBO-A-1",
        tmp_path,
        fig_width=700,
        fig_height=300,
    )

    if output_format == "interactive":
        report = tmp_path / "interactive_metagene_profiling.html"
        assert report.is_file()
        assert "recording-figure" in report.read_text()
        assert figure.image_writes == []
    else:
        assert figure.image_writes == [
            (tmp_path / f"chr1_fiveprime.{output_format}", 700, 300)
        ]


def test_length_cutoff_is_applied_by_annotation_filtering(tmp_path):
    annotation_path = tmp_path / "annotation.gff"
    annotation_path.write_text(
        "chr1\ttest\tCDS\t101\t160\t.\t+\t0\tID=short\n"
        "chr1\ttest\tCDS\t301\t420\t.\t+\t0\tID=long\n"
    )

    starts, stops = annotation.retrieve_annotation_positions(
        annotation_path,
        read_intervals_dict={},
        total_counts_dict={},
        genome_length_dict={"chr1": 1000},
        filtering_methods=["length"],
        mapping_method="fiveprime",
        rpkm_threshold=0,
        overlap_distance=0,
        positions_out_ORF=10,
        positions_in_ORF=20,
        length_cutoff=100,
    )

    assert starts["+"]["chr1"] == [(300, 302)]
    assert stops["+"]["chr1"] == [(417, 419)]


def test_profile_window_geometry_is_transcript_oriented_and_asymmetric():
    plus = annotation.metagene_window_bounds(10, 30, "+", 3, 7)
    minus = annotation.metagene_window_bounds(10, 30, "-", 3, 7)

    assert plus == {"start": (7, 16), "stop": (24, 33)}
    assert minus == {"start": (24, 33), "stop": (7, 16)}


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    "location,start,stop",
    [
        ("left", 3, 20),
        ("right", 20, 38),
    ],
)
def test_annotation_filter_rejects_windows_beyond_contig_edges(
    strand, location, start, stop, tmp_path
):
    annotation_path = tmp_path / f"{strand}-{location}.gff"
    annotation_path.write_text(
        f"chr1\ttest\tCDS\t{start}\t{stop}\t.\t{strand}\t0\tID=edge\n"
    )

    starts, stops = annotation.retrieve_annotation_positions(
        annotation_path,
        read_intervals_dict={},
        total_counts_dict={},
        genome_length_dict={"chr1": 40},
        filtering_methods=["length"],
        mapping_method="fiveprime",
        rpkm_threshold=0,
        overlap_distance=0,
        positions_out_ORF=3,
        positions_in_ORF=7,
        length_cutoff=1,
    )

    assert starts[strand] == {}
    assert stops[strand] == {}


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    "location,start,stop",
    [
        ("left", 4, 20),
        ("right", 20, 37),
    ],
)
def test_annotation_filter_accepts_windows_on_last_valid_contig_coordinate(
    strand, location, start, stop, tmp_path
):
    annotation_path = tmp_path / f"{strand}-{location}.gff"
    annotation_path.write_text(
        f"chr1\ttest\tCDS\t{start}\t{stop}\t.\t{strand}\t0\tID=edge\n"
    )

    starts, stops = annotation.retrieve_annotation_positions(
        annotation_path,
        read_intervals_dict={},
        total_counts_dict={},
        genome_length_dict={"chr1": 40},
        filtering_methods=["length"],
        mapping_method="fiveprime",
        rpkm_threshold=0,
        overlap_distance=0,
        positions_out_ORF=3,
        positions_in_ORF=7,
        length_cutoff=1,
    )

    beginning = start - 1
    end = stop - 1
    expected_start = (beginning, beginning + 2) if strand == "+" else (end - 2, end)
    expected_stop = (end - 2, end) if strand == "+" else (beginning, beginning + 2)
    assert starts[strand]["chr1"] == [expected_start]
    assert stops[strand]["chr1"] == [expected_stop]


@pytest.mark.parametrize(
    "color_list,expected",
    [
        (["#112233", "#abcdef"], ["#112233", "#abcdef", "#112233"]),
        ([], theme.CATEGORICAL[:3]),
    ],
    ids=["configured-colors-cycle", "theme-fallback"],
)
def test_metagene_color_list_controls_profile_traces(
    color_list, expected, tmp_path
):
    coverage = {
        "chr1": {
            read_length: np.arange(1, 7, dtype=np.intp)
            for read_length in (28, 29, 30)
        }
    }

    figures = metagene_profiling.create_metagene_figures(
        coverage,
        {"chr1": {length: values.copy() for length, values in coverage["chr1"].items()}},
        [28, 29, 30],
        tmp_path,
        "fiveprime",
        "raw",
        2,
        4,
        color_list,
    )

    profile = next(
        figure for name, _, figure in figures if name.endswith("(per read length)")
    )
    assert [trace.line.color for trace in profile.data] == expected


def test_tis_advisor_accepts_and_forwards_length_cutoff(monkeypatch):
    args = tis_advisor.parse_arguments(
        [
            "-b",
            "library.bam",
            "-a",
            "annotation.gff",
            "-g",
            "genome.fa",
            "-o",
            "tis-output",
            "--mapping_methods",
            "fiveprime",
            "--length_cutoff",
            "175",
        ]
    )
    captured = []

    class EmptyReader:
        def __init__(self, *args, **kwargs):
            pass

        def output(self):
            return {}, {}

    class EmptyLengthCounter(EmptyReader):
        def output(self):
            return {}

    def retrieve_positions(*call_args):
        captured.append(call_args)
        empty = {"-": {}, "+": {}}
        return empty, empty

    monkeypatch.setattr(tis_advisor.io, "parse_genome_lengths", lambda path: {})
    monkeypatch.setattr(tis_advisor.io, "parse_read_lengths", lambda spec: [25])
    monkeypatch.setattr(tis_advisor, "IntervalReader", EmptyReader)
    monkeypatch.setattr(tis_advisor, "LengthCounter", EmptyLengthCounter)
    monkeypatch.setattr(
        tis_advisor.ann, "retrieve_annotation_positions", retrieve_positions
    )
    monkeypatch.setattr(
        tis_advisor.mg, "metagene_mapping_start", lambda *call_args: {}
    )
    monkeypatch.setattr(
        tis_advisor.mg, "metagene_mapping_stop", lambda *call_args: {}
    )
    monkeypatch.setattr(
        tis_advisor.misc,
        "equalize_dictionary_keys",
        lambda start, stop, *call_args: (start, stop),
    )

    tis_advisor.build_profiles(args)

    assert args.length_cutoff == 175
    assert len(captured) == 1
    assert captured[0][-1] == 175


def test_metagene_rules_quote_settings_and_run_scripts_through_python():
    assert RULES.count("python3 {params.script:q}") == 3
    assert RULES.count("> {log:q} 2>&1") == 3
    assert RULES.count("--length_cutoff") == 2
    assert "--filtering_methods" in RULES
    assert re.search(r"--filtering_method(?!s)", RULES) is None
    assert "{params.colorArgs:q}" in RULES
    assert "if [" not in RULES
    assert "${{colorList}}" not in RULES


@pytest.fixture(scope="module")
def snakemake_command():
    executable = shutil.which("snakemake")
    if executable:
        return [executable]

    conda = shutil.which("conda")
    if conda:
        command = [conda, "run", "-n", "snakemake", "snakemake"]
        probe = subprocess.run(command + ["--version"], capture_output=True, text=True)
        if probe.returncode == 0:
            return command

    pytest.skip("Snakemake is not available for the metagene rule dry-run")


def test_metagene_settings_render_safely_in_a_dry_run(
    snakemake_command,
    genome_file,
    annotation_file,
    samples,
    tmp_path,
):
    workflow_config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    sample_path = tmp_path / "samples.tsv"
    samples[samples["method"] == "RIBO"].iloc[[0]].fillna("").to_csv(
        sample_path, sep="\t", index=False
    )
    workflow_config["biologySettings"].update(
        {
            "genome": str(genome_file),
            "annotation": str(annotation_file),
            "samples": str(sample_path),
        }
    )
    workflow_config["workflowSettings"]["stages"] = ["metagene", "tis_advisor"]
    workflow_config["metageneSettings"].update(
        {
            "lengthCutoff": 175,
            "outputFormats": list(OUTPUT_FORMATS),
            "colorList": ["#123456", "rgb(1, 2, 3)"],
        }
    )

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(workflow_config, sort_keys=False))
    result = subprocess.run(
        [
            *snakemake_command,
            "--dry-run",
            "--printshellcmds",
            "--cores",
            "1",
            "--snakefile",
            str(REPO / "workflow" / "Snakefile"),
            "--directory",
            str(tmp_path),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / ".cache")},
    )
    rendered = result.stdout + result.stderr

    assert result.returncode == 0, rendered
    assert str(REPO / "workflow" / "scripts" / "tis_advisor.py") in rendered
    assert str(REPO / "workflow" / "scripts" / "metagene_profiling.py") in rendered
    assert rendered.count("--length_cutoff 175") == 2
    assert "--filtering_methods overlap length rpkm" in rendered
    assert re.search(r"--filtering_method(?!s)", rendered) is None
    assert "--output_formats interactive svg pdf png jpg" in rendered
    assert "#123456" in rendered
    assert "rgb(1, 2, 3)" in rendered
