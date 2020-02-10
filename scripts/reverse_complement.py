#!/usr/bin/env python
import argparse
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def fasta_reverse_complement(args):
    seqioparse=SeqIO.parse(args.input_fasta_filepath, "fasta")
    outseqio=[]
    for seq_record in seqioparse:
        new_record = SeqRecord(seq_record.seq.reverse_complement(), seq_record.id, seq_record.name, seq_record.description + " reverse complement")
        outseqio.append(new_record)
    return(outseqio)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts fasta to reverse complement.')
    parser.add_argument("--output_fasta_filepath", help='Path to write \
                        fasta output')
    parser.add_argument("--input_fasta_filepath", help='Path to write \
                        fasta input')
    args = parser.parse_args()
    revcomp_records = fasta_reverse_complement(args)
    # write output to fasta file
    SeqIO.write(revcomp_records, args.output_fasta_filepath, "fasta")

if __name__ == '__main__':
    main()
