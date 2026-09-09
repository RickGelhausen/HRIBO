#!/usr/bin/env python
import argparse
import os
from fractions import Fraction

import pysam


def alignment_weight(read, bamfile):
    """Return the fractional contribution of one mapped alignment."""
    if not read.has_tag("NH"):
        return Fraction(1, 1)

    nh = read.get_tag("NH")
    try:
        nh_as_integer = int(nh)
    except (OverflowError, TypeError, ValueError) as error:
        raise ValueError(
            f"{bamfile}: read {read.query_name!r} has invalid NH tag {nh!r}; "
            "expected a positive integer"
        ) from error
    if nh_as_integer <= 0 or nh_as_integer != nh:
        raise ValueError(
            f"{bamfile}: read {read.query_name!r} has invalid NH tag {nh!r}; "
            "expected a positive integer"
        )
    return Fraction(1, nh_as_integer)


def format_mapped_count(value):
    """Keep legacy integer output while representing fractional NH weights."""
    if value.denominator == 1:
        return str(value.numerator)
    return repr(float(value))

def count_mapped_reads(args):
    """Write fractional mapped-read counts and weighted mean read lengths."""

    output_string_mapped = ""
    output_string_length = ""
    bamfiles = sorted(args.bamfiles, key=lambda s: s.lower())
    # run over all input files
    for bamfile in bamfiles:
        total_mapped = {}
        total_length = {}
        with pysam.AlignmentFile(bamfile) as alignment_file:
            for read in alignment_file.fetch():
                if read.is_unmapped:
                    continue

                reference_name = read.reference_name
                weight = alignment_weight(read, bamfile)
                read_length = read.query_length
                if read_length is None:
                    raise ValueError(
                        f"{bamfile}: mapped read {read.query_name!r} has no query length"
                    )
                total_mapped[reference_name] = (
                    total_mapped.get(reference_name, Fraction()) + weight
                )
                total_length[reference_name] = (
                    total_length.get(reference_name, Fraction())
                    + (read_length * weight)
                )

        wildcard = os.path.splitext(os.path.basename(bamfile))[0]

        for key, val in total_mapped.items():
            output_string_mapped += "%s\t%s\t%s\n" % (
                wildcard, key, format_mapped_count(val)
            )

        for key, val in total_length.items():
            output_string_length += "%s\t%s\t%s\n" % (
                wildcard, key, str(float(val / total_mapped[key]))
            )

    with open(args.out_mapped, "w") as f:
        f.write(output_string_mapped)

    with open(args.out_length, "w") as f:
        f.write(output_string_length)

def main():
    parser = argparse.ArgumentParser(
        description="Summarize mapped reads for each BAM and reference sequence."
    )
    parser.add_argument(
        "-b", "--bam", nargs="+", dest="bamfiles", required=True,
        help="Input BAM or SAM files.",
    )
    parser.add_argument(
        "-m", "--out_mapped", dest="out_mapped", required=True,
        help="Output file for mapped-read counts.",
    )
    parser.add_argument(
        "-l", "--out_length", dest="out_length", required=True,
        help="Output file for weighted average read lengths.",
    )
    args = parser.parse_args()

    count_mapped_reads(args)

if __name__ == '__main__':
    main()
