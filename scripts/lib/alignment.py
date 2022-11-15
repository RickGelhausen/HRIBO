import sys

import pysam
import interlap

class IntervalReader():
    """
    Read sam/bam file and create an interlap dictionary for all reads or a selected amount of read lengths
    read_lengths : list (e.g. [30,32,40])
    """
    def __init__(self, alignment_file_path):
        self.alignment_file_path = alignment_file_path

        self.no_accepted_reads_dict = {}
        self.reads_interlap_dict = {}

        print(f"Reading alignment file: {alignment_file_path.stem}")

        self._read_alignment_file()

    def _read_alignment_file(self):
        """
        read the alignment file using pysam
        """

        alignment_file = pysam.AlignmentFile(self.alignment_file_path)
        tmp_dict = {}
        try:
            for read in alignment_file.fetch():
                chrom = read.reference_name

                if read.get_tag("NH") > 1 or read.mapping_quality < 0 or read.is_unmapped:
                    continue

                # start, stop = read.reference_start, read.reference_end-1
                # read_length = stop - start + 1 # alignment length
                #read_length = read.query_length # query read length

                start = read.reference_start
                read_length = read.query_length # query read length
                stop = start + read_length - 1

                strand = "-" if read.is_reverse else "+"

                if chrom in self.no_accepted_reads_dict:
                    self.no_accepted_reads_dict[chrom] += 1
                else:
                    self.no_accepted_reads_dict[chrom] = 1

                interval = (start, stop, read_length)

                if (chrom, strand) in tmp_dict:
                    tmp_dict[(chrom, strand)].append(interval)
                else:
                    tmp_dict[(chrom, strand)] = [interval]

            for key, val in tmp_dict.items():
                inter = interlap.InterLap()
                inter.update(val)
                self.reads_interlap_dict[key] = inter

        except ValueError:
            sys.exit("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.")

    def output(self):
        return self.reads_interlap_dict, self.no_accepted_reads_dict