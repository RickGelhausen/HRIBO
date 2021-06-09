#!/usr/bin/env python
import argparse
import lib.library as hribo


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='HRIBOMerge_offsets')
    parser.add_argument("--in_metagene_directorypath", help='Input metagene directory path', required=True)
    parser.add_argument("--out_filepath", help='Directory path to write output files, if not present the directory will be created', required=True)
    args = parser.parse_args()
    hribo.merge_offset(args.in_metagene_directorypath, args.out_filepath)


if __name__ == '__main__':
    main()
