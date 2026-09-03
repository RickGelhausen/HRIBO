#!/usr/bin/env python3
"""Download and materialize external artifacts only after integrity checks."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path, PurePosixPath

CHUNK_SIZE = 1024 * 1024
CHECKSUM_HEX_LENGTHS = {"md5": 32, "sha256": 64}
INTEGRITY_ERROR_EXIT = 42


class IntegrityError(RuntimeError):
    """Raised when downloaded bytes do not match their recorded identity."""


def _validate_checksum(value: str, algorithm: str) -> str:
    if algorithm not in CHECKSUM_HEX_LENGTHS:
        raise ValueError(f"unsupported checksum algorithm: {algorithm}")
    normalized = value.lower()
    expected_length = CHECKSUM_HEX_LENGTHS[algorithm]
    label = "SHA-256" if algorithm == "sha256" else "MD5"
    if len(normalized) != expected_length:
        raise ValueError(
            f"{label} must contain exactly {expected_length} hexadecimal characters"
        )
    try:
        bytes.fromhex(normalized)
    except ValueError as error:
        raise ValueError(f"{label} must contain only hexadecimal characters") from error
    return normalized


def _new_digest(algorithm: str):
    if algorithm == "md5":
        # UniProt's official HTTPS-served release manifest publishes MD5 for
        # its immutable archive. This is an identity check, not password hashing.
        return hashlib.md5(usedforsecurity=False)
    return hashlib.sha256()


def _temporary_path(destination: Path, suffix: str) -> Path:
    descriptor, name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=suffix,
    )
    os.close(descriptor)
    return Path(name)


def fetch_artifact(
    url: str,
    destination: Path,
    checksum: str,
    *,
    algorithm: str = "sha256",
    expected_size: int | None = None,
    decompress_gzip: bool = False,
    tar_member: str | None = None,
) -> None:
    """Fetch ``url`` and atomically replace ``destination`` after verification.

    The checksum and optional size describe the bytes received from ``url``.
    ``tar_member`` optionally selects one regular file from a verified tar
    archive by its full path or basename. When ``decompress_gzip`` is true, the
    verified download (or selected member) is decompressed into the final
    output. A failed transfer never replaces an existing destination.
    """

    expected_checksum = _validate_checksum(checksum, algorithm)
    if expected_size is not None and expected_size < 0:
        raise ValueError("expected size must be non-negative")

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    downloaded = _temporary_path(destination, ".download")
    materialized: Path | None = None

    try:
        digest = _new_digest(algorithm)
        size = 0
        request = urllib.request.Request(  # noqa: S310 -- caller supplies the URL
            url,
            headers={"User-Agent": "HRIBO verified artifact fetcher"},
        )
        with (
            urllib.request.urlopen(request, timeout=60) as source,  # noqa: S310
            downloaded.open("wb") as target,
        ):
            while chunk := source.read(CHUNK_SIZE):
                target.write(chunk)
                digest.update(chunk)
                size += len(chunk)

        if expected_size is not None and size != expected_size:
            raise IntegrityError(
                f"downloaded size mismatch for {url}: expected {expected_size} bytes, "
                f"received {size}"
            )

        actual_checksum = digest.hexdigest()
        if actual_checksum != expected_checksum:
            label = "SHA-256" if algorithm == "sha256" else "MD5"
            raise IntegrityError(
                f"{label} mismatch for {url}: expected {expected_checksum}, "
                f"received {actual_checksum}"
            )

        if tar_member is not None:
            materialized = _temporary_path(destination, ".materialized")
            with tarfile.open(downloaded, "r:*") as archive:
                matching = [
                    member
                    for member in archive.getmembers()
                    if member.isfile()
                    and (
                        member.name == tar_member
                        or PurePosixPath(member.name).name == tar_member
                    )
                ]
                if len(matching) != 1:
                    raise IntegrityError(
                        f"expected one regular tar member named {tar_member!r}, "
                        f"found {len(matching)}"
                    )
                member_handle = archive.extractfile(matching[0])
                if member_handle is None:
                    raise IntegrityError(f"could not read tar member {matching[0].name!r}")
                with member_handle, materialized.open("wb") as target:
                    if decompress_gzip:
                        with gzip.GzipFile(fileobj=member_handle, mode="rb") as source:
                            shutil.copyfileobj(source, target, length=CHUNK_SIZE)
                    else:
                        shutil.copyfileobj(member_handle, target, length=CHUNK_SIZE)
            os.replace(materialized, destination)
            materialized = None
        elif decompress_gzip:
            materialized = _temporary_path(destination, ".materialized")
            with (
                gzip.open(downloaded, "rb") as source,
                materialized.open("wb") as target,
            ):
                shutil.copyfileobj(source, target, length=CHUNK_SIZE)
            os.replace(materialized, destination)
            materialized = None
        else:
            os.replace(downloaded, destination)
            downloaded = None
    finally:
        if downloaded is not None:
            downloaded.unlink(missing_ok=True)
        if materialized is not None:
            materialized.unlink(missing_ok=True)


def fetch_verified(
    url: str,
    destination: Path,
    sha256: str,
    *,
    expected_size: int | None = None,
    decompress_gzip: bool = False,
) -> None:
    """Compatibility wrapper for the usual SHA-256 download path."""

    fetch_artifact(
        url,
        destination,
        sha256,
        expected_size=expected_size,
        decompress_gzip=decompress_gzip,
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download an external artifact atomically and reject bytes that do not "
            "match the recorded checksum."
        )
    )
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True, type=Path)
    checksum = parser.add_mutually_exclusive_group(required=True)
    checksum.add_argument("--sha256")
    checksum.add_argument("--md5")
    parser.add_argument("--size", type=int)
    parser.add_argument(
        "--tar-member",
        help="materialize this regular file from a verified tar archive",
    )
    parser.add_argument(
        "--decompress-gzip",
        action="store_true",
        help="decompress verified gzip bytes into the final output",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    algorithm = "sha256" if args.sha256 is not None else "md5"
    checksum = args.sha256 if args.sha256 is not None else args.md5
    try:
        fetch_artifact(
            args.url,
            args.output,
            checksum,
            algorithm=algorithm,
            expected_size=args.size,
            decompress_gzip=args.decompress_gzip,
            tar_member=args.tar_member,
        )
    except IntegrityError as error:
        print(f"integrity error: {error}", file=sys.stderr)
        return INTEGRITY_ERROR_EXIT
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
