"""Tests for bugs found while auditing the remaining scripts.

Each test names the behaviour that was wrong, so that a regression is obvious
from the failure rather than from the diff.
"""

import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pandas")

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "workflow" / "scripts"


def run_script(name, *args):
    return subprocess.run(
        [sys.executable, str(SCRIPTS / name), *args],
        capture_output=True, text=True, cwd=str(SCRIPTS),
    )


# --------------------------------------------------------------------------
# samples_to_xlsx
# --------------------------------------------------------------------------


@pytest.fixture
def sample_sheet(tmp_path):
    """Deliberately mixes path shapes, which is what broke the old version."""
    path = tmp_path / "samples.tsv"
    path.write_text(
        "method\tcondition\treplicate\tfastqFile\tfastqFile2\n"
        "RIBO\tA\t1\tfastq/RIBO-A-1_R1.fastq.gz\tfastq/RIBO-A-1_R2.fastq.gz\n"
        "RIBO\tB\t1\t/data/project/fastq/RIBO-B-1_R1.fastq.gz\t\n"
        "RNA\tA\t1\tRNA-A-1_R1.fastq.gz\t\n"
    )
    return path


def read_samples_sheet(path):
    pd = pytest.importorskip("pandas")
    pytest.importorskip("openpyxl")
    return pd.read_excel(path, sheet_name="samples", engine="openpyxl")


def test_samples_to_xlsx_shows_file_names_not_path_fragments(sample_sheet, tmp_path):
    """str[1] took the second path component, so an absolute path became "data"."""
    output = tmp_path / "samples.xlsx"
    result = run_script("samples_to_xlsx.py", "-i", str(sample_sheet), "-o", str(output))
    assert result.returncode == 0, result.stderr

    frame = read_samples_sheet(output)
    assert list(frame["fastqFile"]) == [
        "RIBO-A-1_R1.fastq.gz",
        "RIBO-B-1_R1.fastq.gz",
        "RNA-A-1_R1.fastq.gz",
    ]


def test_samples_to_xlsx_shortens_the_second_read_file(sample_sheet, tmp_path):
    """The column was tested as "Fastqfile2", which never matched "fastqFile2"."""
    output = tmp_path / "samples.xlsx"
    run_script("samples_to_xlsx.py", "-i", str(sample_sheet), "-o", str(output))

    frame = read_samples_sheet(output)
    present = [value for value in frame["fastqFile2"] if isinstance(value, str)]
    assert present == ["RIBO-A-1_R2.fastq.gz"]


def test_samples_to_xlsx_keeps_every_row(sample_sheet, tmp_path):
    """A bare file name used to become NaN rather than the name itself."""
    output = tmp_path / "samples.xlsx"
    run_script("samples_to_xlsx.py", "-i", str(sample_sheet), "-o", str(output))

    frame = read_samples_sheet(output)
    assert len(frame) == 3
    assert frame["fastqFile"].notna().all()


# --------------------------------------------------------------------------
# motif_to_gff
# --------------------------------------------------------------------------


@pytest.fixture
def genomes(tmp_path):
    """A small genome and its reverse complement, as the workflow builds them."""
    Bio = pytest.importorskip("Bio")
    from Bio.Seq import Seq

    # Contains ATG on both strands, so both code paths are exercised.
    sequence = "ATGAAACATTTTATGCAT"
    forward = tmp_path / "genome.fa"
    forward.write_text(f">chr1 test\n{sequence}\n")
    reverse = tmp_path / "genome.rev.fa"
    reverse.write_text(f">chr1 test\n{Seq(sequence).reverse_complement()}\n")
    return forward, reverse, sequence


def run_motif(genomes, tmp_path, motif):
    forward, reverse, _ = genomes
    output = tmp_path / "motifs.gff"
    result = run_script(
        "motif_to_gff.py",
        "--input_genome_fasta_filepath", str(forward),
        "--input_reverse_genome_fasta_filepath", str(reverse),
        "--motif_string", motif,
        "--output_gff3_filepath", str(output),
    )
    assert result.returncode == 0, result.stderr
    rows = [
        line.split("\t")
        for line in output.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    return rows


def test_motif_identifier_format_is_the_same_on_both_strands(genomes, tmp_path):
    """The forward strand emitted a stray colon: ID=chr1:1-3:+: rather than :+."""
    rows = run_motif(genomes, tmp_path, "ATG")
    for row in rows:
        identifier = row[8].split(";")[0]
        assert identifier.endswith(("+", "-")), f"malformed ID: {identifier}"
        assert "::" not in identifier and not identifier.endswith(":")


def test_motif_coordinates_are_correct_on_both_strands(genomes, tmp_path):
    """A hit must actually be the motif at the reported coordinates."""
    from Bio.Seq import Seq

    _, _, sequence = genomes
    for row in run_motif(genomes, tmp_path, "ATG"):
        start, stop, strand = int(row[3]), int(row[4]), row[6]
        segment = sequence[start - 1 : stop]
        if strand == "-":
            segment = str(Seq(segment).reverse_complement())
        assert segment == "ATG", f"{strand} strand {start}-{stop} is {segment!r}"


def test_motif_finds_hits_on_both_strands(genomes, tmp_path):
    strands = {row[6] for row in run_motif(genomes, tmp_path, "ATG")}
    assert strands == {"+", "-"}


# --------------------------------------------------------------------------
# create_reparation_gff
# --------------------------------------------------------------------------


def test_reparation_gff_extends_coordinates_on_both_strands(tmp_path):
    """The strand test used "is" against a literal, which is not guaranteed.

    Both branches extend the ORF by three nucleotides to take in the stop codon,
    outward from the start: downstream on the plus strand, upstream on the minus.
    """
    predicted = tmp_path / "Predicted_ORFs.txt"
    columns = [
        "ORF_locus", "strand", "length", "start_codon", "ribo_count", "ribo_rpkm",
        "ribo_coverage", "SD_score", "SD_pos", "prob", "ORF_type", "Reference",
        "Distance_from_aTIS",
    ]
    rows = [
        ["chr1:100-199", "+", "100", "ATG", "10", "1.0", "0.9", "5", "-8", "0.9", "sORF", "aTIS", "0"],
        ["chr1:300-399", "-", "100", "ATG", "10", "1.0", "0.9", "5", "-8", "0.9", "sORF", "aTIS", "0"],
    ]
    predicted.write_text(
        "\t".join(columns) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n"
    )

    output = tmp_path / "out.gff"
    result = run_script(
        "create_reparation_gff.py", "-c", "A", "-r", "1",
        "-i", str(predicted), "-o", str(output),
    )
    assert result.returncode == 0, result.stderr

    coordinates = {}
    for line in output.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        coordinates[fields[6]] = (int(fields[3]), int(fields[4]))

    assert coordinates["+"] == (100, 202), "plus strand stop codon not added"
    assert coordinates["-"] == (297, 399), "minus strand stop codon not added"


def test_reparation_gff_emits_no_syntax_warning():
    """`is` against a string literal is a SyntaxWarning, and not guaranteed."""
    source = (SCRIPTS / "create_reparation_gff.py").read_text()
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error", SyntaxWarning)
        compile(source, "create_reparation_gff.py", "exec")
