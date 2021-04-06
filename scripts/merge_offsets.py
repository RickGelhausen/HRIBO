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
        profiling_type_dir_key=os.path.basename(os.path.normpath(profiling_type_dir))
        offset_dict[profiling_type_dir_key]={}
        norm_type_dirs = [f.path for f in os.scandir(profiling_type_dir) if f.is_dir()]
        for norm_type_dir in norm_type_dirs:
            #samples
            norm_type_dir_key=os.path.basename(os.path.normpath(norm_type_dir))
            offset_dict[profiling_type_dir_key][norm_type_dir_key]={}
            samples = [f.path for f in os.scandir(norm_type_dir) if f.is_dir()]
            for sample in samples:
                sample_key=os.path.basename(os.path.normpath(sample))
                offset_dict[profiling_type_dir_key][norm_type_dir_key][sample_key]={}
                path=sample #os.path.dirname(sample)
                files = [f for f in listdir(sample) if isfile(join(sample, f))]
                json_files = [f for f in files if f.find(".json") != -1]
                for json_file in json_files:
                    json_filename = os.path.basename(json_file)
                    json_filename = re.sub('\_length_offset.json$', '', json_filename)
                    offset_dict[profiling_type_dir_key][norm_type_dir_key][sample_key][json_filename]={}
                    #print(json_filename)
                    with open((path + "/" +  json_file), 'r') as jf:
                        data=jf.read()
                        read_dict = json.loads(data)
                        offset_dict[profiling_type_dir_key][norm_type_dir_key][sample_key][json_filename]=read_dict
    #delete readlengthstat node obtained from adding length_reads_dict
    offset_dict.pop('readlengthstats', None)
    with open((out_path + "/" + 'merged_offsets.json'), 'w') as fp:
        json.dump(offset_dict, fp)

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
