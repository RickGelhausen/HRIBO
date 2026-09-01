#!/usr/bin/env python
import argparse
import re

from Bio import SeqIO


def parse_motifs(motif_string):
    """Return non-empty, literal DNA motifs in a consistent case."""
    return [motif.strip().upper() for motif in motif_string.split(",") if motif.strip()]


def motif_spans(sequence, motif):
    """Find every literal occurrence, including overlapping occurrences."""
    return [
        (match.start(), match.start() + len(motif))
        for match in re.finditer(f"(?={re.escape(motif)})", sequence.upper())
    ]

#HRIBO/scripts/motif_to_gff.py --input_genome_fasta_filepath genomes/genome.fa --input_reverse_genome_fasta_filepath genomes/genome.revfa --motif_string ATG --output_gff3_filepath test.gff3
def motif_gff3_forward_strand(args):
    seqioparse = SeqIO.parse(args.input_genome_fasta_filepath, "fasta")
    out_entries = ""
    counter = 1
    motifs = parse_motifs(args.motif_string)
    for seq_record in seqioparse:
        for motif in motifs:
            coordinates = motif_spans(str(seq_record.seq), motif)
            for (start_coordinate, end_coordinate) in coordinates:
                frame = start_coordinate % 3
                start_out = start_coordinate + 1
                end_out = end_coordinate
                id_out = f"{seq_record.id}:{start_out}-{end_out}:+"
                motif_entry = f"{seq_record.id}\tHRIBO\tnucleotide_motif\t{start_out}\t{end_out}"\
                            + f"\t.\t+\t.\tID={id_out};Name={motif};frame={frame};\n"

                out_entries += motif_entry
                counter += 1
    return out_entries

def motif_gff3_reverse_strand(args):
    seqioparse = SeqIO.parse(args.input_reverse_genome_fasta_filepath, "fasta")
    out_entries = ""
    counter = 1
    motifs = parse_motifs(args.motif_string)
    for seq_record in seqioparse:
        for motif in motifs:
            reverse_sequence = str(seq_record.seq)
            length_reverse_seq = len(reverse_sequence)
            coordinates = motif_spans(reverse_sequence, motif)
            for (start_coordinate, end_coordinate) in coordinates:
                frame = (length_reverse_seq - end_coordinate) % 3
                start_out = length_reverse_seq - end_coordinate + 1
                end_out = length_reverse_seq - start_coordinate
                id_out = f"{seq_record.id}:{start_out}-{end_out}:-"
                motif_entry = f"{seq_record.id}\tHRIBO\tnucleotide_motif\t{start_out}\t{end_out}"\
                            + f"\t.\t-\t.\tID={id_out};Name={motif};frame={frame};\n"

                out_entries += motif_entry
                counter += 1
    return out_entries


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Searches the genome for motifs and builds gff3 tracks for hits.')
    parser.add_argument(
        "--motif_string",
        nargs="?",
        const="",
        default="",
        help="Comma-separated motif strings; an empty value produces a header-only track",
    )
    parser.add_argument("--output_gff3_filepath", help='Path to write gff3 output')
    parser.add_argument("--input_genome_fasta_filepath", help='Path to read genome fasta input')
    parser.add_argument("--input_reverse_genome_fasta_filepath", help='Path to read reverse genome fasta input')
    args = parser.parse_args()
    plusgff3 = motif_gff3_forward_strand(args)
    minusgff3 = motif_gff3_reverse_strand(args)
    gff3 = str(plusgff3+minusgff3)

    with open(args.output_gff3_filepath, "w") as f:
        f.write("##gff-version 3\n")
        f.write(gff3)

if __name__ == '__main__':
    main()
