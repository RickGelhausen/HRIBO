#!/usr/bin/env python
import argparse
import os
import pandas as pd
import json
from os import listdir
from os.path import isfile, join
import re

def merge_offset(in_metagene_directorypath,out_path):
    #tis/tts
    offset_dict = {}
    profiling_type_dirs = [f.path for f in os.scandir(in_metagene_directorypath) if f.is_dir()]
    for profiling_type_dir in profiling_type_dirs:
        #norm/raw
        norm_type_dirs = [f.path for f in os.scandir(profiling_type_dir) if f.is_dir()]
        for norm_type_dir in norm_type_dirs:
            #samples
            samples = [f.path for f in os.scandir(norm_type_dir) if f.is_dir()]
            for sample in samples:
                files = [f for f in listdir(sample) if isfile(join(sample, f))]
                json_files = [f for f in files if f.find(".json") != -1]
                for json_file in json_files:
                    json_filename = os.path.basename(json_file)
                    json_filename = re.sub('\.json$', '', json_filename)
                    print(json_filename)

    return ""


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='HRIBOMerge_offsets')
    parser.add_argument("--in_metagene_directorypath", help='Input metagene directory path', required=True)
    parser.add_argument("--out_filepath", help='Directory path to write output files, if not present the directory will be created', required=True)
    args = parser.parse_args()
    merge_offset(args.in_metagene_directorypath, args.out_filepath)


if __name__ == '__main__':
    main()
