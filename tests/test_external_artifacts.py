"""Integrity and reproducibility tests for downloaded workflow artifacts."""

import gzip
import hashlib
import io
import math
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

from fetch_verified import IntegrityError, fetch_artifact, fetch_verified


REPO = Path(__file__).resolve().parent.parent
DEEPRIBO_RULES = (REPO / "workflow" / "rules" / "deepribo.smk").read_text()
DELTATE_RULES = (REPO / "workflow" / "rules" / "diffex_deltate.smk").read_text()
REPARATION_RULES = (REPO / "workflow" / "rules" / "reparation.smk").read_text()
FETCHER = REPO / "workflow" / "scripts" / "fetch_verified.py"


def deepribo_parameter_reader():
    """Load the ordinary Python preamble before Snakemake's rule syntax."""

    preamble = DEEPRIBO_RULES.split("rule deepriboGetModel:", 1)[0]
    namespace = {"math": math}
    exec(compile(preamble, "deepribo.smk", "exec"), namespace)
    return namespace["read_parameters"]


def sha256(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()


def test_deepribo_parameter_reader_accepts_one_finite_bounded_pair(tmp_path):
    parameters = tmp_path / "parameters.txt"
    parameters.write_text(" 12.5 , 0.75 \n")
    reader = deepribo_parameter_reader()

    assert reader(parameters, 0) == "12.5"
    assert reader(parameters, 1) == "0.75"


def test_deepribo_parameter_reader_rejects_malformed_cutoffs(tmp_path):
    reader = deepribo_parameter_reader()
    invalid = (
        "",
        "1",
        "1,0.5,extra",
        "NA,0.5",
        "NaN,0.5",
        "Inf,0.5",
        "-1,0.5",
        "1,-0.1",
        "1,1.1",
        "1,0.5\n2,0.6",
    )

    for index, content in enumerate(invalid):
        parameters = tmp_path / f"invalid-{index}.txt"
        parameters.write_text(content)
        with pytest.raises(ValueError):
            reader(parameters, 0)


def test_fetch_verified_materializes_exact_bytes(tmp_path):
    contents = b"stable model bytes\x00\x01"
    source = tmp_path / "source model.pt"
    destination = tmp_path / "download directory" / "model.pt"
    source.write_bytes(contents)

    fetch_verified(
        source.as_uri(),
        destination,
        sha256(contents),
        expected_size=len(contents),
    )

    assert destination.read_bytes() == contents
    assert list(destination.parent.glob(f".{destination.name}.*")) == []


def test_fetch_verified_decompresses_only_after_verification(tmp_path):
    contents = b">protein one\nMPEPTIDE\n"
    source = tmp_path / "swiss prot.fasta.gz"
    destination = tmp_path / "database" / "uniprot_sprot.fasta"
    with gzip.GzipFile(filename=source, mode="wb", mtime=0) as handle:
        handle.write(contents)
    compressed = source.read_bytes()

    fetch_verified(
        source.as_uri(),
        destination,
        sha256(compressed),
        expected_size=len(compressed),
        decompress_gzip=True,
    )

    assert destination.read_bytes() == contents
    assert list(destination.parent.glob(f".{destination.name}.*")) == []


def test_verified_tar_member_can_supply_an_archived_fallback(tmp_path):
    contents = b">archived protein\nMPEPTIDE\n"
    compressed = gzip.compress(contents, mtime=0)
    source = tmp_path / "uniprot release.tar.gz"
    with tarfile.open(source, "w:gz") as archive:
        member = tarfile.TarInfo("release/files/uniprot_sprot.fasta.gz")
        member.size = len(compressed)
        archive.addfile(member, io.BytesIO(compressed))
    archive_bytes = source.read_bytes()
    archive_md5 = hashlib.md5(archive_bytes, usedforsecurity=False).hexdigest()
    destination = tmp_path / "database" / "uniprot_sprot.fasta"

    fetch_artifact(
        source.as_uri(),
        destination,
        archive_md5,
        algorithm="md5",
        expected_size=len(archive_bytes),
        tar_member="uniprot_sprot.fasta.gz",
        decompress_gzip=True,
    )

    assert destination.read_bytes() == contents
    assert list(destination.parent.glob(f".{destination.name}.*")) == []


@pytest.mark.parametrize(
    ("expected_hash", "expected_size", "message"),
    [
        ("0" * 64, None, "SHA-256 mismatch"),
        (sha256(b"new bytes"), 999, "downloaded size mismatch"),
    ],
)
def test_integrity_failure_preserves_previous_output(
    tmp_path, expected_hash, expected_size, message
):
    source = tmp_path / "source"
    destination = tmp_path / "output"
    source.write_bytes(b"new bytes")
    destination.write_bytes(b"previous verified bytes")

    with pytest.raises(IntegrityError, match=message):
        fetch_verified(
            source.as_uri(),
            destination,
            expected_hash,
            expected_size=expected_size,
        )

    assert destination.read_bytes() == b"previous verified bytes"
    assert list(tmp_path.glob(f".{destination.name}.*")) == []


def test_corrupt_gzip_preserves_previous_output(tmp_path):
    source = tmp_path / "broken.gz"
    destination = tmp_path / "database.fasta"
    compressed = b"\x1f\x8bnot a complete gzip stream"
    source.write_bytes(compressed)
    destination.write_bytes(b">previous\nMPEPTIDE\n")

    with pytest.raises((EOFError, gzip.BadGzipFile)):
        fetch_verified(
            source.as_uri(),
            destination,
            sha256(compressed),
            decompress_gzip=True,
        )

    assert destination.read_bytes() == b">previous\nMPEPTIDE\n"
    assert list(tmp_path.glob(f".{destination.name}.*")) == []


def test_cli_distinguishes_integrity_mismatch_from_transfer_failure(tmp_path):
    source = tmp_path / "model.pt"
    destination = tmp_path / "downloaded.pt"
    source.write_bytes(b"unexpected bytes")

    result = subprocess.run(
        [
            sys.executable,
            str(FETCHER),
            "--url",
            source.as_uri(),
            "--sha256",
            "0" * 64,
            "--output",
            str(destination),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 42
    assert "integrity error: SHA-256 mismatch" in result.stderr
    assert not destination.exists()


def test_external_artifacts_are_content_addressed():
    combined = DEEPRIBO_RULES + DELTATE_RULES + REPARATION_RULES

    assert "docker://gelhausr/deepribo:latest" not in combined
    assert "docker://gelhausr/deltate:latest" not in combined
    assert "DeepRibo/raw/master" not in combined
    assert "@sha256:" in DEEPRIBO_RULES
    assert "@sha256:" in DELTATE_RULES
    assert "quay.io/biocontainers/reparation_blast@sha256:" in REPARATION_RULES
    assert (
        "6852b3b69b532039a5d674115b9cbfd2953f6479cf187bcc769e2d899fcbc288"
        in REPARATION_RULES
    )
    assert "../envs/reparation.yaml" not in REPARATION_RULES
    assert not (REPO / "workflow/envs/reparation.yaml").exists()
    assert "fetch_verified.py" in DEEPRIBO_RULES
    assert "fetch_verified.py" in REPARATION_RULES
    assert "--sha256" in DEEPRIBO_RULES
    assert "--sha256" in REPARATION_RULES
    assert "release 2026_02" in REPARATION_RULES
    assert "previous_releases/" in REPARATION_RULES
    assert "release-2026_02/" in REPARATION_RULES
    assert "--md5" in REPARATION_RULES
    assert "--tar-member" in REPARATION_RULES
    assert "compact_status" in REPARATION_RULES
