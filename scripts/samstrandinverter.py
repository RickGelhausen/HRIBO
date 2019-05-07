#!/usr/bin/env python
import csv
import sys
import argparse

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Swaps strand in sam file.')
    parser.add_argument("--sam_in_filepath", default=None, help='Path to read sam input')
    parser.add_argument("--sam_out_filepath", default=None, help='Path to write sam output')
    args = parser.parse_args()
    insamhandle = open(args.sam_in_filepath, 'r')
    insam = csv.reader(insamhandle, dialect="excel-tab")
    outsamhandle = open(args.sam_out_filepath, 'w')
    outsam = csv.writer(outsamhandle, dialect="excel-tab")
    for entry in insam :
        #Don't try to modify the header!!!
        if(entry[0][0] == "@") :
            outsam.writerow(entry)
            continue
        #Swap strand
        if((int(entry[1]) & 0x10) == 0x10) :
            entry[1] = int(entry[1]) - 0x10
        else :
            entry[1] = int(entry[1]) + 0x10
        outsam.writerow(entry)

if __name__ == '__main__':
    main()
