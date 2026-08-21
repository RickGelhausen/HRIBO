#!/usr/bin/env python
"""
Write per-strand A-site occupancy bedgraph files for DeepRibo.

DeepRibo maps each read to a single genomic position using "a 12 nt offset from
the 3' end of the read". The 3' end is the rightmost aligned base on the forward
strand and the leftmost aligned base on the reverse strand, and the offset moves
into the read in both cases.

Coordinates here are 0-based throughout, matching the bedgraph format that
DataParser.py consumes alongside the bedtools genomecov output.

Author: Rick Gelhausen
"""

import argparse
import sys
from collections import Counter

import pysam

# Distance from the 3' end of a read to the ribosomal A-site, per DeepRibo.
A_SITE_OFFSET = 12


class ASiteOccupancy:
    """Count A-site positions per strand from an alignment file."""

    def __init__(self, alignment_file, prefix, offset=A_SITE_OFFSET):
        self.alignment_file = alignment_file
        self.offset = offset
        self.forward_out = f"{prefix}_asite_fwd.bedgraph"
        self.reverse_out = f"{prefix}_asite_rev.bedgraph"

        self.forward = Counter()
        self.reverse = Counter()
        self.counted = 0
        self.skipped = 0
        self.out_of_bounds = 0

    def a_site_position(self, reference_start, read_length, is_reverse):
        """The A-site of one read, as a 0-based genomic coordinate.

        forward: the 3' end is the rightmost base, so the offset moves left
        reverse: the 3' end is the leftmost base, so the offset moves right
        """
        if is_reverse:
            return reference_start + self.offset
        return reference_start + read_length - 1 - self.offset

    def run(self):
        with pysam.AlignmentFile(self.alignment_file) as alignments:
            lengths = dict(zip(alignments.references, alignments.lengths))
            for read in alignments.fetch():
                # Flags are a bit field: testing `flag == 0` or `flag == 16`, as
                # the previous version did, silently drops every read that also
                # carries the secondary, supplementary or duplicate bit.
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                    self.skipped += 1
                    continue

                read_length = read.query_length or read.infer_query_length()
                if not read_length:
                    self.skipped += 1
                    continue

                position = self.a_site_position(
                    read.reference_start, read_length, read.is_reverse
                )

                contig_length = lengths.get(read.reference_name)
                if position < 0 or (contig_length is not None and position >= contig_length):
                    self.out_of_bounds += 1
                    continue

                target = self.reverse if read.is_reverse else self.forward
                target[(read.reference_name, position)] += 1
                self.counted += 1

        self._write(self.forward, self.forward_out)
        self._write(self.reverse, self.reverse_out)

    @staticmethod
    def _write(counts, path):
        """Write a bedgraph, sorted by contig and position.

        Sorted output matters because the file sits beside the bedtools
        genomecov coverage tracks, which are sorted, and because unsorted
        bedgraph is invalid for most consumers.
        """
        with open(path, "w") as handle:
            for (contig, position), value in sorted(counts.items()):
                handle.write(f"{contig}\t{position}\t{position + 1}\t{value}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Write A-site occupancy bedgraph files for DeepRibo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--alignment_file", required=True,
                        help="Alignment file to read (sam/bam).")
    parser.add_argument("--output_file_prefix", dest="output_prefix", required=True,
                        help="Prefix for the two output bedgraph files.")
    parser.add_argument("--offset", type=int, default=A_SITE_OFFSET,
                        help="Distance from the 3' end of a read to the A-site.")
    args = parser.parse_args()

    occupancy = ASiteOccupancy(args.alignment_file, args.output_prefix, args.offset)
    occupancy.run()

    print(
        f"A-site positions written for {occupancy.counted} reads "
        f"({len(occupancy.forward)} forward, {len(occupancy.reverse)} reverse positions).",
        file=sys.stderr,
    )
    if occupancy.skipped:
        print(f"Skipped {occupancy.skipped} unmapped or non-primary alignments.", file=sys.stderr)
    if occupancy.out_of_bounds:
        print(
            f"Discarded {occupancy.out_of_bounds} A-site positions falling outside their contig.",
            file=sys.stderr,
        )
    if occupancy.counted == 0:
        sys.exit(
            "No A-site positions were derived from "
            f"{args.alignment_file}. DeepRibo cannot run on an empty A-site track."
        )


if __name__ == "__main__":
    main()
