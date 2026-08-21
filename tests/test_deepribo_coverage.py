"""Tests for the DeepRibo A-site occupancy track.

DeepRibo maps a read to a single position using a 12 nt offset from its 3' end.
The previous implementation computed the reverse-strand position from
`read.pos - read_length`, which places the A-site roughly a full read length
outside the alignment, so every reverse-strand gene fed DeepRibo a misplaced
signal. These tests pin the geometry down on both strands.
"""

import subprocess
import sys
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "workflow" / "scripts" / "coverage_deepribo.py"

CONTIG = "chr1"
CONTIG_LENGTH = 5000
OFFSET = 12


def write_bam(path, reads):
    """reads: iterable of (reference_start, length, is_reverse[, flag_extra])."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}],
    }
    records = sorted(reads, key=lambda r: r[0])
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for index, entry in enumerate(records):
            start, length, is_reverse = entry[:3]
            extra_flag = entry[3] if len(entry) > 3 else 0
            segment = pysam.AlignedSegment()
            segment.query_name = f"r{index}"
            segment.query_sequence = "A" * length
            segment.query_qualities = pysam.qualitystring_to_array("I" * length)
            segment.flag = (16 if is_reverse else 0) | extra_flag
            segment.reference_id = 0
            segment.reference_start = start
            segment.mapping_quality = 255
            segment.cigartuples = [(0, length)]
            segment.set_tag("NH", 1)
            out.write(segment)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    Path(unsorted).unlink()


def run_script(bam, prefix):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--alignment_file", str(bam),
         "--output_file_prefix", str(prefix)],
        capture_output=True, text=True,
    )
    return result


def read_bedgraph(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        contig, start, end, value = line.split("\t")
        rows.append((contig, int(start), int(end), int(value)))
    return rows


@pytest.fixture
def prefix(tmp_path):
    return tmp_path / "asite"


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------


@pytest.mark.parametrize("start,length", [(1000, 30), (2000, 28), (500, 35), (10, 24)])
def test_forward_a_site_is_offset_from_the_three_prime_end(tmp_path, prefix, start, length):
    bam = tmp_path / "fwd.bam"
    write_bam(bam, [(start, length, False)])
    run_script(bam, prefix)

    rows = read_bedgraph(f"{prefix}_asite_fwd.bedgraph")
    three_prime_end = start + length - 1  # rightmost base on the forward strand
    assert rows == [(CONTIG, three_prime_end - OFFSET, three_prime_end - OFFSET + 1, 1)]


@pytest.mark.parametrize("start,length", [(1000, 30), (2000, 28), (500, 35), (10, 24)])
def test_reverse_a_site_is_offset_from_the_three_prime_end(tmp_path, prefix, start, length):
    bam = tmp_path / "rev.bam"
    write_bam(bam, [(start, length, True)])
    run_script(bam, prefix)

    rows = read_bedgraph(f"{prefix}_asite_rev.bedgraph")
    three_prime_end = start  # leftmost base on the reverse strand
    assert rows == [(CONTIG, three_prime_end + OFFSET, three_prime_end + OFFSET + 1, 1)]


@pytest.mark.parametrize("is_reverse", [False, True])
@pytest.mark.parametrize("length", [24, 28, 30, 35])
def test_a_site_always_falls_inside_the_read(tmp_path, prefix, is_reverse, length):
    """The regression this file exists for: the A-site left the alignment."""
    start = 1000
    bam = tmp_path / "read.bam"
    write_bam(bam, [(start, length, is_reverse)])
    run_script(bam, prefix)

    suffix = "rev" if is_reverse else "fwd"
    rows = read_bedgraph(f"{prefix}_asite_{suffix}.bedgraph")
    assert len(rows) == 1
    position = rows[0][1]
    assert start <= position <= start + length - 1, (
        f"A-site {position} lies outside the read span [{start}, {start + length - 1}]"
    )


def test_strands_are_written_to_separate_files(tmp_path, prefix):
    bam = tmp_path / "both.bam"
    write_bam(bam, [(1000, 30, False), (1000, 30, True)])
    run_script(bam, prefix)

    forward = read_bedgraph(f"{prefix}_asite_fwd.bedgraph")
    reverse = read_bedgraph(f"{prefix}_asite_rev.bedgraph")
    assert len(forward) == 1 and len(reverse) == 1
    assert forward[0][1] == 1000 + 30 - 1 - OFFSET
    assert reverse[0][1] == 1000 + OFFSET


def test_reads_at_the_same_a_site_are_summed(tmp_path, prefix):
    bam = tmp_path / "stack.bam"
    write_bam(bam, [(1000, 30, False)] * 5)
    run_script(bam, prefix)

    rows = read_bedgraph(f"{prefix}_asite_fwd.bedgraph")
    assert rows == [(CONTIG, 1017, 1018, 5)]


# --------------------------------------------------------------------------
# Flag handling
# --------------------------------------------------------------------------


def test_secondary_and_duplicate_alignments_are_skipped(tmp_path, prefix):
    bam = tmp_path / "flags.bam"
    write_bam(bam, [
        (1000, 30, False),            # primary, counted
        (1100, 30, False, 0x100),     # secondary
        (1200, 30, False, 0x400),     # duplicate
        (1300, 30, False, 0x800),     # supplementary
    ])
    run_script(bam, prefix)

    rows = read_bedgraph(f"{prefix}_asite_fwd.bedgraph")
    assert rows == [(CONTIG, 1017, 1018, 1)]


# --------------------------------------------------------------------------
# Output form
# --------------------------------------------------------------------------


def test_output_is_sorted_by_position(tmp_path, prefix):
    """Read lengths vary, so A-site order does not follow alignment order."""
    bam = tmp_path / "mixed.bam"
    write_bam(bam, [(1000, 35, False), (1005, 24, False), (1002, 30, False)])
    run_script(bam, prefix)

    positions = [row[1] for row in read_bedgraph(f"{prefix}_asite_fwd.bedgraph")]
    assert positions == sorted(positions)


def test_bedgraph_intervals_are_single_bases(tmp_path, prefix):
    bam = tmp_path / "single.bam"
    write_bam(bam, [(1000, 30, False), (2000, 28, True)])
    run_script(bam, prefix)

    for suffix in ("fwd", "rev"):
        for _, start, end, _ in read_bedgraph(f"{prefix}_asite_{suffix}.bedgraph"):
            assert end == start + 1


def test_out_of_bounds_a_sites_are_discarded(tmp_path, prefix):
    """A read at the very start of a contig can push a reverse A-site past its end."""
    bam = tmp_path / "edge.bam"
    write_bam(bam, [(CONTIG_LENGTH - 5, 4, True)])
    result = run_script(bam, prefix)

    rows = read_bedgraph(f"{prefix}_asite_rev.bedgraph")
    assert rows == []
    assert "outside their contig" in result.stderr


def test_empty_result_fails_loudly(tmp_path, prefix):
    """An empty A-site track would make DeepRibo fail much later and obscurely."""
    bam = tmp_path / "empty.bam"
    write_bam(bam, [(1000, 30, False, 0x100)])
    result = run_script(bam, prefix)

    assert result.returncode != 0
    assert "empty A-site track" in result.stderr


# --------------------------------------------------------------------------
# Configurable offset
#
# DeepRibo's published 12 nt was derived from E. coli, so the offset has to be
# adjustable for other organisms and digestion protocols.
# --------------------------------------------------------------------------


def run_script_with_offset(bam, prefix, offset):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--alignment_file", str(bam),
         "--output_file_prefix", str(prefix), "--offset", str(offset)],
        capture_output=True, text=True,
    )


@pytest.mark.parametrize("offset", [0, 6, 12, 15, 18])
def test_offset_shifts_both_strands_symmetrically(tmp_path, prefix, offset):
    start, length = 1000, 30
    bam = tmp_path / "both.bam"
    write_bam(bam, [(start, length, False), (start, length, True)])
    run_script_with_offset(bam, prefix, offset)

    forward = read_bedgraph(f"{prefix}_asite_fwd.bedgraph")
    reverse = read_bedgraph(f"{prefix}_asite_rev.bedgraph")

    assert forward[0][1] == start + length - 1 - offset
    assert reverse[0][1] == start + offset


def test_default_offset_is_twelve(tmp_path, prefix):
    """The DeepRibo published value stays the default."""
    bam = tmp_path / "default.bam"
    write_bam(bam, [(1000, 30, False)])
    run_script(bam, prefix)
    explicit_prefix = tmp_path / "explicit"
    run_script_with_offset(bam, explicit_prefix, 12)

    assert read_bedgraph(f"{prefix}_asite_fwd.bedgraph") == \
           read_bedgraph(f"{explicit_prefix}_asite_fwd.bedgraph")


def test_zero_offset_lands_on_the_three_prime_end(tmp_path, prefix):
    start, length = 1000, 30
    bam = tmp_path / "zero.bam"
    write_bam(bam, [(start, length, False), (start, length, True)])
    run_script_with_offset(bam, prefix, 0)

    assert read_bedgraph(f"{prefix}_asite_fwd.bedgraph")[0][1] == start + length - 1
    assert read_bedgraph(f"{prefix}_asite_rev.bedgraph")[0][1] == start
