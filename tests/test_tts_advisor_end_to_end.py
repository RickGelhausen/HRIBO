"""Termination advice must recover a stop-bound A-site from genomic read positions.

The planted peaks use the first stop-codon base on each strand and contain no
initiation peak. These tests cover the complete BAM-to-report path, so coordinate
or library-type mistakes cannot pass merely by giving a scoring helper its own
preferred input coordinates.
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pysam")
pytest.importorskip("interlap")

import simulate_riboseq as sim

REPO = Path(__file__).resolve().parent.parent
ADVISOR = REPO / "workflow" / "scripts" / "tis_advisor.py"


def run_advisor(reference_dir, bam, out_dir, *, library_type=None,
                read_lengths="24-35", filtering_methods=("length",), rpkm_threshold=0):
    command = [
        sys.executable, str(ADVISOR),
        "-b", str(bam),
        "-a", str(reference_dir / "annotation.gff"),
        "-g", str(reference_dir / "genome.fa"),
        "-o", str(out_dir),
        "-r", read_lengths,
        "--rpkm_threshold", str(rpkm_threshold),
        "--filtering_methods", *filtering_methods,
        "--include_plotly_js", "online",
    ]
    if library_type is not None:
        command.extend(["--library_type", library_type])
    result = subprocess.run(command, capture_output=True, text=True,
                            cwd=str(REPO / "workflow" / "scripts"))
    assert result.returncode == 0, result.stderr
    return json.loads((out_dir / "tis_recommendation.json").read_text())


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    directory = tmp_path_factory.mktemp("tts-reference")
    sim.write_genome(directory / "genome.fa")
    sim.write_annotation(directory / "annotation.gff")
    return directory


@pytest.mark.parametrize("periodic", [True, False])
@pytest.mark.parametrize(
    "read_end,expected_offset",
    [("fiveprime", sim.PLANTED_TTS_A_SITE_OFFSET),
     ("threeprime", sim.PLANTED_TTS_THREE_PRIME_A_SITE_OFFSET)],
)
def test_recovers_stop_a_site_and_protocol_read_end(
    reference, tmp_path, periodic, read_end, expected_offset,
):
    bam = tmp_path / "TTS-condition-1.bam"
    sim.write_bam(bam, periodic=periodic, seed=11, peak="stop", anchor=read_end)
    payload = run_advisor(reference, bam, tmp_path / "out")

    recommendation = payload["recommendation"]
    assert payload["library_type"] == "TTS"
    assert payload["anchor"] == recommendation["anchor"] == "stop"
    assert payload["site"] == recommendation["site"] == "A"
    assert payload["chosen_read_end"] == recommendation["read_end"] == read_end
    assert recommendation["read_lengths"], "the stop-only library has a clear peak"
    assert set(recommendation["read_lengths"]) <= sim.GOOD_LENGTHS
    assert set(recommendation["offsets"].values()) == {expected_offset}
    assert recommendation["confidence"] in {"high", "medium"}
    assert set(payload["read_ends"]) == {"fiveprime", "threeprime"}
    for comparison in payload["read_ends"].values():
        assert comparison["recommendation"]["anchor"] == "stop"
        assert comparison["recommendation"]["site"] == "A"
        assert all(score["anchor"] == "stop" and score["site"] == "A"
                   for score in comparison["read_lengths"])


@pytest.mark.parametrize("strand", ["+", "-"])
def test_first_stop_base_is_calibrated_on_each_strand(reference, tmp_path, strand):
    """Pooling two strands must not conceal an incorrectly oriented stop anchor."""
    bam = tmp_path / "TTS-condition-1.bam"
    sim.write_bam(bam, periodic=False, seed=17, peak="stop", strands=(strand,))
    payload = run_advisor(reference, bam, tmp_path / "out")

    recommendation = payload["recommendation"]
    assert recommendation["read_lengths"]
    assert recommendation["read_end"] == "fiveprime"
    assert set(recommendation["offsets"].values()) == {sim.PLANTED_TTS_A_SITE_OFFSET}


def test_strong_initiation_signal_cannot_drive_tts_advice(reference, tmp_path):
    bam = tmp_path / "TTS-condition-1.bam"
    sim.write_bam(bam, periodic=False, seed=2, peak="start")
    payload = run_advisor(reference, bam, tmp_path / "out")
    initiation = run_advisor(reference, bam, tmp_path / "start-control",
                             library_type="TIS")

    assert initiation["recommendation"]["read_lengths"]
    assert set(initiation["recommendation"]["offsets"].values()) == {sim.PLANTED_OFFSET}
    assert payload["recommendation"]["read_lengths"] == []
    assert payload["recommendation"]["confidence"] == "none"
    assert payload["recommendation"]["warnings"]
    assert payload["chosen_read_end"] is None


def test_explicit_tts_type_overrides_an_arbitrary_filename(reference, tmp_path):
    bam = tmp_path / "my_termination_experiment.bam"
    sim.write_bam(bam, periodic=False, seed=11, peak="stop")
    payload = run_advisor(reference, bam, tmp_path / "out", library_type="TTS")

    assert payload["library_type"] == "TTS"
    assert payload["recommendation"]["read_lengths"]
    assert set(payload["recommendation"]["offsets"].values()) == {
        sim.PLANTED_TTS_A_SITE_OFFSET
    }


@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_reports_label_measured_a_site_and_derived_p_site(reference, tmp_path, read_end):
    bam = tmp_path / "TTS-condition-1.bam"
    sim.write_bam(bam, periodic=False, seed=11, peak="stop", anchor=read_end)
    out = tmp_path / "out"
    payload = run_advisor(reference, bam, out)

    # The output filenames remain stable for workflow/report consumers.
    assert (out / "tis_recommendation.json").is_file()
    html = (out / "tis_recommendation.html").read_text()
    assert "TTS peak advice" in html
    assert ">A-site offset (nt)</th>" in html
    assert ">Derived P-site offset (nt)</th>" in html
    assert ">Derived A-site offset (nt)</th>" not in html
    assert "A-site" in html
    assert "stop codon" in html.lower()
    assert "psiteOffsets:" not in html
    assert "Suggested configuration" not in html
    assert "DeepRibo A-site offset advice" not in html
    assert payload["deepribo_a_site"]["suggested_offset"] is None
    assert payload["deepribo_a_site"]["applied"] is False
    recommendation = payload["recommendation"]
    assert recommendation["a_site_offsets"] == recommendation["offsets"]
    assert (out / "orfbounder/manifest.json").is_file()
    assert "orfbounder/" in html
    for end, comparison in payload["read_ends"].items():
        end_recommendation = comparison["recommendation"]
        directory = out / "orfbounder" / end
        if not end_recommendation["read_lengths"]:
            assert not (directory / "offsets.json").exists()
            continue
        assert json.loads((directory / "read_lengths.json").read_text()) == {
            bam.stem: ",".join(str(length) for length in sorted(end_recommendation["read_lengths"]))
        }
        direction = -1 if end == "fiveprime" else 1
        # Termination exports use the measured A-site without deriving a P-site.
        assert json.loads((directory / "offsets.json").read_text()) == {
            bam.stem: {
                length: direction * offset
                for length, offset in end_recommendation["offsets"].items()
            }
        }

    with (out / "read_length_evidence.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert {"offset", "anchor", "site", "read_end", "a_site_offset",
                "p_site_offset"} <= set(reader.fieldnames)
        rows = {int(row["read_length"]): row for row in reader}
    assert set(rows) == set(sim.ALL_LENGTHS)
    for length in recommendation["read_lengths"]:
        row = rows[length]
        assert row["anchor"] == "stop"
        assert row["site"] == "A"
        assert row["read_end"] == read_end
        a_offset = int(row["a_site_offset"])
        assert int(row["offset"]) == a_offset
        # A and P differ by one codon along the transcript. The distance from
        # a 5' end decreases, while the distance back from a 3' end increases.
        expected_p = a_offset - 3 if read_end == "fiveprime" else a_offset + 3
        assert int(row["p_site_offset"]) == expected_p
        assert recommendation["p_site_offsets"][str(length)] == expected_p


def test_long_termination_footprints_are_not_limited_to_monomer_offsets(
    reference, tmp_path, monkeypatch,
):
    """Configured longer footprints still expose their stop-bound A-site."""
    monkeypatch.setattr(sim, "ALL_LENGTHS", [58, 59, 60])
    monkeypatch.setattr(sim, "GOOD_LENGTHS", {58, 59, 60})
    bam = tmp_path / "TTS-disome-1.bam"
    sim.write_bam(bam, periodic=False, seed=11, peak="stop", anchor="threeprime",
                  three_prime_offset=13)
    payload = run_advisor(reference, bam, tmp_path / "out", read_lengths="58-60")

    recommendation = payload["recommendation"]
    assert recommendation["read_lengths"]
    assert recommendation["read_end"] == "threeprime"
    assert set(recommendation["offsets"].values()) == {13}
    fiveprime_offsets = {
        score["read_length"]: score["offset"]
        for score in payload["read_ends"]["fiveprime"]["read_lengths"]
        if score["usable"]
    }
    assert fiveprime_offsets
    assert all(offset == length - 1 - 13
               for length, offset in fiveprime_offsets.items())


@pytest.mark.parametrize("anchor,site,axis", [("start", "P", "x"), ("stop", "A", "x2")])
@pytest.mark.parametrize("read_end,direction", [("fiveprime", -1), ("threeprime", 1)])
def test_offset_markers_use_the_calibrated_boundary_panel(anchor, site, axis,
                                                         read_end, direction):
    """Both panels remain visible, but offsets belong to their calibrated codon."""
    import numpy as np
    import pandas as pd
    from lib import plotting

    coordinates = np.arange(-25, 26)
    frame = pd.DataFrame({"coordinates": coordinates, "28": np.ones(coordinates.size)})
    figure = plotting.plot_metagene_heatmap(
        frame, frame, [28], "Metagene", offsets={28: 15}, read_end=read_end,
        anchor=anchor, site=site,
    )

    markers = [trace for trace in figure.data if trace.type == "scatter"]
    assert len(markers) == 1
    assert markers[0].xaxis == axis
    assert list(markers[0].x) == [direction * 15]
    assert markers[0].name == f"estimated {site}-site offset"


@pytest.mark.parametrize("read_end", ["fiveprime", "threeprime"])
def test_tts_plots_align_the_a_site_with_the_first_stop_base(
    reference, tmp_path, monkeypatch, read_end,
):
    """Inspect figures from real BAM profiling, including pooled peak geometry."""
    import numpy as np
    import tis_advisor as advisor

    bam = tmp_path / "TTS-condition-1.bam"
    sim.write_bam(bam, periodic=False, seed=11, peak="stop", anchor=read_end)
    out = tmp_path / "out"
    figures = {}
    render_report = advisor.render_report

    def capture_report(*args, **kwargs):
        figures.update(dict(args[3]))
        return render_report(*args, **kwargs)

    monkeypatch.setattr(advisor, "render_report", capture_report)
    monkeypatch.setattr(sys, "argv", [
        str(ADVISOR), "-b", str(bam),
        "-a", str(reference / "annotation.gff"),
        "-g", str(reference / "genome.fa"),
        "-o", str(out), "-r", "24-35",
        "--filtering_methods", "length", "--include_plotly_js", "online",
    ])
    advisor.main()
    payload = json.loads((out / "tis_recommendation.json").read_text())

    heatmap = figures["Read length against position"]
    markers = next(trace for trace in heatmap.data if trace.type == "scatter")
    assert markers.xaxis == "x2"
    assert markers.name == "estimated A-site offset"
    direction = -1 if read_end == "fiveprime" else 1
    for length, offset in zip(markers.y, markers.customdata):
        stop_panel = heatmap.data[1]
        row = list(stop_panel.y).index(length)
        peak = np.argmax(np.asarray(stop_panel.z)[row])
        assert stop_panel.x[peak] == direction * offset

    pooled = figures["Pooled, offset corrected"]
    trace = pooled.data[0]
    assert trace.x[np.argmax(trace.y)] == 0
    assert trace.name == "pooled A-sites"
    assert "stop" in trace.hovertemplate and "A-sites" in trace.hovertemplate
    assert pooled.layout.xaxis.title.text.startswith("Distance from stop codon")
    assert pooled.layout.yaxis.title.text == "Pooled A-sites"
    assert any(annotation.text == "stop codon" for annotation in pooled.layout.annotations)
    assert payload["recommendation"]["site"] == "A"


def test_rpkm_filter_retains_stop_footprints_with_three_prime_ends_outside_cds(
    reference, tmp_path,
):
    """Termination footprints overlap CDS even when their measured end does not."""
    bam = tmp_path / "TTS-condition-1.bam"
    sim.write_bam(bam, periodic=False, seed=11, peak="stop", anchor="threeprime",
                  body_count=0, noise=False)
    payload = run_advisor(reference, bam, tmp_path / "out",
                          filtering_methods=("rpkm", "length"), rpkm_threshold=10)

    recommendation = payload["recommendation"]
    assert recommendation["read_lengths"], "overlapping stop footprints failed the RPKM filter"
    assert recommendation["read_end"] == "threeprime"
    assert set(recommendation["offsets"].values()) == {
        sim.PLANTED_TTS_THREE_PRIME_A_SITE_OFFSET
    }
