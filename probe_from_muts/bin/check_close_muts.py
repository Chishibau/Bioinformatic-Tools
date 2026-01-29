#!/usr/bin/env python

"""
Created on Mon Jan 06 2026

@author: zshong
"""

import os
import pandas as pd
import argparse
import sys
import numpy as np

def close_muts(path_output, path_input, size):
    if os.path.isfile(path_input):
        muts=pd.read_csv(path_input, sep="\t")
    else:
        print("Error: please provide valid file path to muts with 'Position' column")
        return
    
    out_df = pd.DataFrame()

    # get muts at exact same location
    duplicate_rows_mask = muts.duplicated(subset=['Position'], keep=False)
    out_df = muts[duplicate_rows_mask]

    # get muts that are close to each other
    no_dup = muts.drop_duplicates(subset=['Position'], keep=False)

    # set tolerance for proximity as half the probe length
    tolerance = size / 2

    # create arrays to compare each mut position
    positions = no_dup['Position'].values
    x, y = np.meshgrid(positions, positions)

    close_pos_mask = (np.abs(x - y) <= tolerance) & (x != y)
    close_pairs = set()
    for i, j in zip(*np.where(close_pos_mask)):
        # Store as a sorted tuple to avoid duplicates
        pair = tuple(sorted((positions[i], positions[j])))
        close_pairs.add(pair)

    # convert close mut pairs into individual rows
    close_positions = [pos for pair in close_pairs for pos in pair]
    close_df = no_dup[no_dup['Position'].isin(close_positions)]
    out_df = pd.concat([out_df, close_df], ignore_index=True)

    # sort by position
    out_df = out_df.sort_values('Position').reset_index(drop=True)

    # get muts that are NOT close to others
    result_df = pd.merge(muts, out_df, how='left', indicator=True)
    result_df = result_df[result_df['_merge'] == 'left_only']
    result_df = result_df.drop(columns=['_merge']).reset_index(drop=True)

    # save output files
    path_filtered = path_output + "/no_close_muts.tsv"
    path_close_muts = path_output + "/close_loc_muts.tsv"

    result_df.to_csv(path_filtered, index=False, sep='\t')
    out_df.to_csv(path_close_muts, index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser(description="get muts that are within probe-length of each other")
    parser.add_argument("--path_output", required=True, help="abs path to output file")
    parser.add_argument("--path_input", required=True, help="abs path to input file")
    parser.add_argument("--size", required=False, type=int, default=120, help="length of probes")
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()

    close_muts(args.path_output, args.path_input, args.size)
    
if __name__ == "__main__":
    main()