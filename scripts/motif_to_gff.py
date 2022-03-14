#!/usr/bin/env python
import argparse
import re

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#HRIBO/scripts/motif_to_gff.py --input_genome_fasta_filepath genomes/genome.fa --input_reverse_genome_fasta_filepath genomes/genome.revfa --motif_string ATG --output_gff3_filepath test.gff3
def motif_gff3_forward_strand(args):
    seqioparse = SeqIO.parse(args.input_genome_fasta_filepath, "fasta")
    out_entries = ""
    counter = 1
    motifs = args.motif_string.split(",")
    for seq_record in seqioparse:
        for motif in motifs:
            coordinates = [m.span() for m in re.finditer(str(motif), str(seq_record.seq))]
            for (start_coordinate, end_coordinate) in coordinates:
                frame = start_coordinate % 3
                start_out = start_coordinate + 1
                end_out = end_coordinate
                id_out = f"{seq_record.id}:{start_out}-{end_out}:+"
                motif_entry = f"{seq_record.id}\tHRIBO\tnucleotide_motif\t{start_out}\t{end_out}"\
                            + f"\t.\t+\t.\tID={id_out}:;Name={motif};Frame={frame};\n"

                out_entries += motif_entry
                counter += 1
    return out_entries

def motif_gff3_reverse_strand(args):
    seqioparse = SeqIO.parse(args.input_reverse_genome_fasta_filepath, "fasta")
    out_entries = ""
    counter = 1
    motifs = args.motif_string.split(",")
    for seq_record in seqioparse:
        for motif in motifs:
            reverse_sequence = str(seq_record.seq)
            length_reverse_seq = len(reverse_sequence)
            coordinates = [m.span() for m in re.finditer(str(motif), reverse_sequence)]
            for (start_coordinate, end_coordinate) in coordinates:
                frame = (length_reverse_seq - end_coordinate) % 3
                start_out = length_reverse_seq - end_coordinate + 1
                end_out = length_reverse_seq - start_coordinate
                id_out = f"{seq_record.id}:{start_out}-{end_out}:-"
                motif_entry = f"{seq_record.id}\tHRIBO\tnucleotide_motif\t{start_out}\t{end_out}"\
                            + f"\t.\t-\t.\tID={id_out};Name={motif};Frame={frame};\n"

                out_entries += motif_entry
                counter += 1
    return out_entries


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Searches the genome for motifs and builds gff3 tracks for hits.')
    parser.add_argument("--motif_string", help='Comma separated motif strings')
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
