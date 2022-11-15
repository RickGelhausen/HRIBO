#!/usr/bin/env python
# samtools view -b -F 4 data/TIS-TolC-2.bam > data/TIS-TolC-2_mapped.bam
# samtools index data/TIS-TolC-2_mapped.bam
# samtools faidx exp5genome.fa
# time ./MetageneProfiling.py --in_bam_filepath ecoli/TIS-TolC-2.bam --in_gff_filepath ecoli/annotation.gff --out_plot_filepath ecoli/metagene
import argparse
import lib.library as hribo


def main():

    # store commandline args
    parser = argparse.ArgumentParser(description='HRIBOMetageneProfiling')
    parser.add_argument("--in_bam_filepath", help='Input bam path', required=True)
    parser.add_argument("--in_gff_filepath", help='Input gff path', required=True)
    parser.add_argument("--in_readlengthstat_filepath", help='Input read length statistics json path', required=False)
    parser.add_argument("--cpu_cores", help='Number of cpu cores to use', type=int, default=1)
    parser.add_argument("--min_read_length", help='Minimal read length to consider', type=int, default=27)
    parser.add_argument("--max_read_length", help='Maximal read length to consider', type=int, default=33)
    parser.add_argument("--out_plot_filepath", help='Directory path to write output files, if not present the directory will be created', required=True)
    parser.add_argument("--normalization", help='Toggles normalization by average read count per nucleotide', action='store_true')
    parser.add_argument("--noise_reduction_analysis", help='Toggles noise reduction analysis', action='store_true')
    parser.add_argument("--input_type", help='Input type, either TIS or TTS.', default="TIS")
    args = parser.parse_args()
    hribo.meta_geneprofiling_p(args.input_type, args.in_gff_filepath, args.in_bam_filepath, args.out_plot_filepath, args.cpu_cores, args.min_read_length, args.max_read_length, args.normalization, args.noise_reduction_analysis, args.in_readlengthstat_filepath)


if __name__ == '__main__':
    main()
