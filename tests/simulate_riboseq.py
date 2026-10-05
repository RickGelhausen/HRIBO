"""Simulate a Ribo-seq BAM with a known ribosomal site, for end-to-end testing.

Reads are placed so that their 5' end sits PLANTED_OFFSET nt upstream of each
annotated start codon (the initiation peak), plus an elongation background along
the ORF body, plus uniform noise. Read lengths in GOOD_LENGTHS get the signal;
the rest get noise only. ``peak="stop"`` instead plants a termination peak with
the A-site on the first base of the annotated stop codon.
"""

import random
import sys
from pathlib import Path

import pysam

PLANTED_OFFSET = 12
# Used when anchoring on the 3' end: the 3' end sits this far downstream of the
# start codon, so the 5' end is the one that spreads out by read length.
PLANTED_THREE_PRIME_OFFSET = 15
# A terminating ribosome has its A-site one codon downstream of its P-site.
PLANTED_TTS_A_SITE_OFFSET = PLANTED_OFFSET + 3
PLANTED_TTS_THREE_PRIME_A_SITE_OFFSET = PLANTED_THREE_PRIME_OFFSET - 3
GOOD_LENGTHS = {28, 29, 30}
ALL_LENGTHS = list(range(24, 36))
CONTIGS = {"NC_000913.3": 60000, "pPlasmid1": 20000}


def genes(periodic=True):
    """Same layout as the validation fixture: CDS every 450 nt, 300 nt long."""
    out = []
    index = 0
    for name, length in CONTIGS.items():
        pos = 200
        while pos + 400 < length:
            index += 1
            feature = "rRNA" if index % 7 == 0 else ("tRNA" if index % 11 == 0 else "CDS")
            strand = "+" if index % 2 else "-"
            out.append((name, feature, pos, pos + 299, strand, index))
            pos += 450
    return out


def write_annotation(path):
    lines = ["##gff-version 3"]
    for name, feature, start, end, strand, index in genes():
        lines.append("\t".join([
            name, "sim", feature, str(start), str(end), ".", strand, "0",
            f"ID=gene{index};locus_tag=b{index:04d}",
        ]))
    Path(path).write_text("\n".join(lines) + "\n")


def write_genome(path):
    rng = random.Random(0)
    with open(path, "w") as fh:
        for name, length in CONTIGS.items():
            fh.write(f">{name} simulated\n")
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            for i in range(0, length, 70):
                fh.write(seq[i:i + 70] + "\n")


def write_bam(path, periodic=True, signal=True, seed=0, anchor="fiveprime",
              three_prime_offset=None, peak="start", five_prime_offset=None,
              strands=("+", "-"), body_count=64, noise=True):
    """Plant a start P-site or stop A-site peak at a defined read end.

    "fiveprime" mimics a protocol where the 5' end is precisely defined, so the 3'
    end spreads out by read length. "threeprime" is the opposite, which is the
    case in many bacteria where nuclease digestion trims the 3' end sharply.
    Stop-codon coordinates are calculated directly from GFF CDS bounds, without
    using the advisor's metagene axis or offset conversion code.
    """
    if peak not in {"start", "stop"}:
        raise ValueError("peak must be 'start' or 'stop'")
    if anchor not in {"fiveprime", "threeprime"}:
        raise ValueError("anchor must be 'fiveprime' or 'threeprime'")
    if five_prime_offset is None:
        five_prime_offset = (
            PLANTED_TTS_A_SITE_OFFSET if peak == "stop" else PLANTED_OFFSET
        )
    if three_prime_offset is None:
        three_prime_offset = (
            PLANTED_TTS_THREE_PRIME_A_SITE_OFFSET
            if peak == "stop" else PLANTED_THREE_PRIME_OFFSET
        )
    rng = random.Random(seed)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in CONTIGS.items()],
    }
    records = []

    for name, feature, start, end, strand, _ in genes():
        if feature != "CDS" or strand not in strands:
            continue
        # GFF is 1-based inclusive; pysam positions are 0-based.
        cds_start = start - 1
        cds_end = end - 1
        if peak == "stop":
            # The stop codon occupies the final three bases in transcript order;
            # its FIRST base is end-2 on '+' and beginning+2 on '-'.
            target = cds_end - 2 if strand == "+" else cds_start + 2
        else:
            target = cds_start if strand == "+" else cds_end

        def read_position(site, read_length, read_end):
            if read_end == "fiveprime":
                if strand == "+":
                    return site - five_prime_offset
                return site + five_prime_offset - read_length + 1
            if strand == "+":
                return site + three_prime_offset - read_length + 1
            return site - three_prime_offset

        for read_length in ALL_LENGTHS:
            carries_signal = signal and read_length in GOOD_LENGTHS

            if carries_signal:
                # Boundary peak, anchored on whichever end the protocol defines.
                for _ in range(rng.randint(14, 22)):
                    pos = read_position(target, read_length, anchor)
                    records.append((name, pos, read_length, strand))

                # elongation signal along the body
                for codon in range(6, 6 + body_count):
                    if periodic and rng.random() > 0.55:
                        continue
                    if not periodic and rng.random() > 0.18:
                        continue
                    body = codon * 3 + (0 if periodic else rng.randint(0, 2))
                    if peak == "stop":
                        body = -body
                    site = target + body if strand == "+" else target - body
                    # Preserve the legacy start fixture's elongation background.
                    body_read_end = anchor if peak == "stop" else "fiveprime"
                    pos = read_position(site, read_length, body_read_end)
                    records.append((name, pos, read_length, strand))

            # uniform noise across the whole locus, for every read length
            for _ in range(rng.randint(2, 6) if noise else 0):
                pos = rng.randint(cds_start - 150, cds_end + 150)
                records.append((name, pos, read_length, strand))

    records = [r for r in records if 0 <= r[1] and r[1] + r[2] < CONTIGS[r[0]]]
    records.sort(key=lambda r: (list(CONTIGS).index(r[0]), r[1]))

    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for i, (name, pos, read_length, strand) in enumerate(records):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * read_length
            a.query_qualities = pysam.qualitystring_to_array("I" * read_length)
            a.flag = 16 if strand == "-" else 0
            a.reference_id = list(CONTIGS).index(name)
            a.reference_start = pos
            a.mapping_quality = 255
            a.cigartuples = [(0, read_length)]
            a.set_tag("NH", 1)
            out.write(a)

    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    Path(unsorted).unlink()
    return len(records)


if __name__ == "__main__":
    outdir = Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)
    write_genome(outdir / "genome.fa")
    write_annotation(outdir / "annotation.gff")
    n1 = write_bam(outdir / "RIBO-A-1.bam", periodic=True, signal=True, seed=1)
    n2 = write_bam(outdir / "RIBO-A-2.bam", periodic=False, signal=True, seed=2)
    n3 = write_bam(outdir / "RNA-A-1.bam", periodic=False, signal=False, seed=3)
    print(f"planted offset {PLANTED_OFFSET}, good lengths {sorted(GOOD_LENGTHS)}")
    print(f"reads: periodic={n1} nonperiodic={n2} nosignal={n3}")
