#!/usr/bin/env python
"""
Created on Fri Oct 17 14:07:07 2025

@author: zshong
"""

import os
import pandas as pd
import argparse
import sys

def generate_bed(path_tsv, path_output, chrom_col, pos_col, id_col, size):
    data = pd.read_csv(path_tsv, sep='\t')
    
    # generate regions of 200 bp around each mutation position (+- 100 bp upstream and downstream)
    bed = pd.DataFrame()
    bed['chrom'] = data[chrom_col]
    bed['start'] = data[pos_col] - size
    bed['end'] = data[pos_col] + size
    bed['name'] = data[id_col] + "_" + data[chrom_col] + "_" + data[pos_col].astype(str)
    
    width = str(size*2)
    path_file = path_output + "/all_muts_" + width + "bp_regions.bed"

    # Remove old file if it exists
    try:
        os.remove(path_file)
    except FileNotFoundError:
        pass

    # generate bed file and save in target dir
    bed.to_csv(path_file, index=False, sep='\t', header=None)

    print("RESULT_PATH: " + path_file)

def main():
    parser = argparse.ArgumentParser(description="generate bed region file from single mut positions")
    parser.add_argument("--path_output", required=True, help="abs path to output dir")
    parser.add_argument("--path_tsv", required=True, help="path to tsv muts table")
    parser.add_argument("--chrom_col", required=False, default='CHROM', help="name of chr column in data sheet")
    parser.add_argument("--pos_col", required=False,  default='POSITION', help="name of mut position column in data sheet")
    parser.add_argument("--id_col", required=False,  default='patient', help="name of id column in data sheet")
    parser.add_argument("--size", required=False, type=int, default=100, help="number of bp to include before and after mut pos")
    args = parser.parse_args()

    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()

    generate_bed(args.path_tsv, args.path_output, args.chrom_col, args.pos_col, args.id_col, args.size)
    
if __name__ == "__main__":
    main()