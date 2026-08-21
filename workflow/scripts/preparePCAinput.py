#!/usr/bin/env python
import pandas as pd
import argparse
import numpy as np

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Create the meta sheet for PCA analysis.")
    parser.add_argument("-s", "--sample_sheet", action="store", dest="sample_sheet", required=True, help="Sample sheet.")
    parser.add_argument("-o", "--output_file", action="store", dest="output_file", required=True, help= "The output folder.")
    args = parser.parse_args()


    samples = pd.read_csv(args.sample_sheet, dtype=str, sep="\t")
    print(samples)
    samples = samples.sort_values(by=["method", "condition", "replicate"], ascending=True, key=lambda x: x if np.issubdtype(x.dtype, np.number) else x.str.lower())
    print(samples)
    out_string = "\tsampletype\texpr\n"
    for line in samples.itertuples():
        out_string += f"{line.method}_{line.condition}_{line.replicate}\t{line.method}_{line.condition}\t{line.condition}\n"

    with open(args.output_file, "w") as f:
        f.write(out_string)

if __name__ == '__main__':
    main()