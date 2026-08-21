#!/usr/bin/env python
import pysam
import re
import sys
import argparse

class generateASiteOccupancy:

    def __init__(self, alignment_file, prefix):
        self.alignment_file = alignment_file # sam / bam format

        self.a_site_s_dict = {}
        self.a_site_as_dict = {}

        self.asite_sense_out = prefix + "_asite_fwd.bedgraph"
        self.asite_antisense_out = prefix + "_asite_rev.bedgraph"

        self._fill_dictionary()

    def _write_output(self):
        """
        writes the contents of dictionary to file
        """
        with open(self.asite_sense_out, "w") as of:
            # go through every entry in the dictionary
            for key, val in self.a_site_s_dict.items():
                # write line to file
                of.write("%s\t%s\t%s\t%s\n" % (key[0], key[1], key[2], val))

        with open(self.asite_antisense_out, "w") as of:
            # go through every entry in the dictionary
            for key, val in self.a_site_as_dict.items():
                # write line to file
                of.write("%s\t%s\t%s\t%s\n" % (key[0], key[1], key[2], val))

    def _a_site_occupancy_fwd(self, flag, region, stop):
        """
        calculate the a-site occupancy for the forward strand
        """
        a_site_start = stop - 12
        a_site_stop = a_site_start + 1

        entry = (region, a_site_start, a_site_stop)
        if entry in self.a_site_s_dict:
            self.a_site_s_dict[entry] += 1
        else:
            self.a_site_s_dict[entry] = 1

    def _a_site_occupancy_rev(self, flag, region, stop):
        """
        calculate the a-site occupancy for the reverse strand
        """
        a_site_start = stop + 11
        a_site_stop = a_site_start + 1

        entry = (region, a_site_start, a_site_stop)
        if entry in self.a_site_as_dict:
            self.a_site_as_dict[entry] += 1
        else:
            self.a_site_as_dict[entry] = 1

    def _fill_dictionary(self):
        """
        main function handling the different possible cases of calculating the coverage and a-site
        """
        # read bam or sam file
        samfile = pysam.AlignmentFile(self.alignment_file)
        for read in samfile.fetch():
            # get required attributes
            flag = int(read.flag) # flag: 0-forward-strand 4-unmapped 16-reverse-strand
            reference_name = read.reference_name # identifier for the region
            reference_pos = read.pos + 1 # 0- based leftmost mapping position (reference_pos)
            read_length = len(read.query_sequence)

            if flag == 0:
                self._a_site_occupancy_fwd(flag, reference_name, reference_pos+read_length)
            elif flag == 16:
                self._a_site_occupancy_rev(flag, reference_name, reference_pos-read_length)
            else:
                continue

        samfile.close()
        self._write_output()

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='get the occupancy of a site as bed files')
    parser.add_argument("--alignment_file", action="store", dest="alignment_file", required=True
                        , help="the alignment file to be used (sam/bam)")
    parser.add_argument("--output_file_prefix", action="store", dest="output_prefix", required=True
                        , help="the prefix for the output file.")
    args = parser.parse_args()

    gBG = generateASiteOccupancy(args.alignment_file, args.output_prefix)
    print("DONE!")

if __name__ == '__main__':
    main()
