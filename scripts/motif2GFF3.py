#!/usr/bin/env python
import argparse
import re
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#SPtools/scripts/motif2GFF3.py --input_genome_fasta_filepath genomes/genome.fa --input_reverse_genome_fasta_filepath genomes/genome.revfa --motif_string ATG --output_gff3_filepath test.gff3
def motif_gff3_plus_strand(args):
  seqioparse=SeqIO.parse(args.input_genome_fasta_filepath, "fasta")
  outentries=""
  for seq_record in seqioparse:
    coordinates = [m.span() for m in re.finditer(str(args.motif_string),str(seq_record.seq))]
    for (start_coordinate,end_coordinate) in coordinates:
        motifentry = seq_record.id + "\t" + "." + "\t" + "nucleotide_motif" + "\t" + str(start_coordinate) + "\t" +  str(end_coordinate)  + "\t" + "." + "\t" + "+" + "\t" + "." + "\t" + "Note=motif;\n"
        outentries+=motifentry
  return(outentries)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Searches the genome for motifs and builds gff3 tracks for hits.')
    parser.add_argument("--motif_string", help='Motif string')
    parser.add_argument("--output_gff3_filepath", help='Path to write gff3 output')
    parser.add_argument("--input_genome_fasta_filepath", help='Path to read genome fasta input')
    parser.add_argument("--input_reverse_genome_fasta_filepath", help='Path to read reverse genome fasta input')
    args = parser.parse_args()
    plusgff3 = motif_gff3_plus_strand(args)
    gff3=str(plusgff3)
    f = open(args.output_gff3_filepath, 'wt', encoding='utf-8')
    f.write(gff3) 
if __name__ == '__main__':
    main()

