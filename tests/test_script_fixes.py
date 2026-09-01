"""Tests for bugs found while auditing the remaining scripts.

Each test names the behaviour that was wrong, so that a regression is obvious
from the failure rather than from the diff.
"""

import shutil
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


def assert_valid_gff3(path):
    genome_tools = shutil.which("gt")
    if genome_tools is None:
        return
    result = subprocess.run(
        [genome_tools, "gff3validator", str(path)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


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
    pytest.importorskip("Bio")
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


def test_motif_input_is_case_insensitive_and_empty_input_has_no_hits(genomes, tmp_path):
    assert run_motif(genomes, tmp_path, "atg") == run_motif(genomes, tmp_path, "ATG")
    assert run_motif(genomes, tmp_path, "") == []


def test_motif_option_without_a_shell_token_produces_a_header_only_track(
    genomes, tmp_path
):
    """Snakemake omits the token when an empty configured list is formatted."""
    forward, reverse, _ = genomes
    output = tmp_path / "no-alternative-starts.gff"

    result = run_script(
        "motif_to_gff.py",
        "--input_genome_fasta_filepath",
        str(forward),
        "--input_reverse_genome_fasta_filepath",
        str(reverse),
        "--motif_string",
        "--output_gff3_filepath",
        str(output),
    )

    assert result.returncode == 0, result.stderr
    assert output.read_text() == "##gff-version 3\n"


def test_motif_search_retains_overlapping_hits(tmp_path):
    pytest.importorskip("Bio")
    from Bio.Seq import Seq

    forward = tmp_path / "overlap.fa"
    reverse = tmp_path / "overlap.rev.fa"
    forward.write_text(">chr1\nAAAA\n")
    reverse.write_text(f">chr1\n{Seq('AAAA').reverse_complement()}\n")

    plus_rows = [row for row in run_motif((forward, reverse, "AAAA"), tmp_path, "AAA") if row[6] == "+"]
    assert [(int(row[3]), int(row[4])) for row in plus_rows] == [(1, 3), (2, 4)]


def test_alternative_start_track_uses_the_configured_codons():
    rules = (
        Path(__file__).resolve().parent.parent
        / "workflow"
        / "rules"
        / "visualization.smk"
    ).read_text()
    block = rules.split("rule alternativeStartCodonTrack:", 1)[1].split(
        "rule stopCodonTrack:", 1
    )[0]

    assert '",".join(str(codon).upper() for codon in CODONS)' in block
    assert "GTG,TTG,CTG" not in block
    assert "{params.motifs:q}" in block


# --------------------------------------------------------------------------
# enrich_annotation
# --------------------------------------------------------------------------


def test_enrichment_reports_a_malformed_id_without_nameerror(tmp_path):
    annotation = tmp_path / "malformed.gff"
    output = tmp_path / "enriched.gff"
    annotation.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t9\t.\t+\t.\tbadid=value;\n"
    )

    result = run_script(
        "enrich_annotation.py",
        "-a",
        str(annotation),
        "-o",
        str(output),
    )

    assert result.returncode != 0
    assert "Missing ID in annotation row" in result.stderr
    assert "NameError" not in result.stderr


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
    phases = []
    for line in output.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        coordinates[fields[6]] = (int(fields[3]), int(fields[4]))
        phases.append(fields[7])

    assert coordinates["+"] == (100, 202), "plus strand stop codon not added"
    assert coordinates["-"] == (297, 399), "minus strand stop codon not added"
    assert phases == ["0", "0"]
    assert output.read_text().startswith("##gff-version 3\n")
    assert "orf_type=sORF" in output.read_text()
    assert "ORF_type=" not in output.read_text()
    assert_valid_gff3(output)


def test_reparation_gff_emits_no_syntax_warning():
    """`is` against a string literal is a SyntaxWarning, and not guaranteed."""
    source = (SCRIPTS / "create_reparation_gff.py").read_text()
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error", SyntaxWarning)
        compile(source, "create_reparation_gff.py", "exec")


def test_deepribo_gff_stores_distance_outside_phase(tmp_path):
    predicted = tmp_path / "predictions.csv"
    columns = [
        "filename",
        "filename_counts",
        "label",
        "in_gene",
        "strand",
        "coverage",
        "coverage_elo",
        "rpk",
        "rpk_elo",
        "start_site",
        "start_codon",
        "stop_site",
        "stop_codon",
        "locus",
        "prot_seq",
        "nuc_seq",
        "pred",
        "pred_rank",
        "SS",
        "dist",
        "SS_pred_rank",
    ]
    values = [
        "sample",
        "counts",
        "True",
        "False",
        "+",
        "1",
        "1",
        "1",
        "1",
        "10",
        "ATG",
        "20",
        "TAA",
        "chr1:10-20",
        "M",
        "ATG",
        "0.9",
        "1",
        "1",
        "-1",
        "1",
    ]
    predicted.write_text(",".join(columns) + "\n" + ",".join(values) + "\n")

    output = tmp_path / "deepribo.gff"
    result = run_script(
        "create_deepribo_gff.py",
        "-c",
        "A",
        "-r",
        "1",
        "-i",
        str(predicted),
        "-o",
        str(output),
    )
    assert result.returncode == 0, result.stderr

    row = [line for line in output.read_text().splitlines() if not line.startswith("#")][0]
    fields = row.split("\t")
    assert fields[7] == "0"
    assert "deepribo_distance=-1" in fields[8]
    assert "condition=A" in fields[8]
    assert "Condition=" not in fields[8]
    assert_valid_gff3(output)
