#!/usr/bin/env python
"""
Created on Mon Oct 20 13:24:40 2025

@author: zshong
"""

import os
import pandas as pd
import argparse
import sys

def closest_center_probes(path_output, path_input, size):
    if os.path.isfile(path_input):
        probes=pd.read_csv(path_input, sep="\t")
    else:
        print("Error: please provide valid file path to txt with candidate probes")
    
    probes.columns = ["probe_id"]

    # split probe_id into individual comps, chr num, mut pos, etc.
    data = probes['probe_id'].str.split('_', expand=True)
    # data[0] = data[0] + '_' + data[1]
    # data[1] = data[2]
    # data[2] = data[3]
    # data.drop(3, axis=1, inplace=True)
    data[3] = data[2].str.split('::').str[1]
    data[3] = data[3].str.split(':').str[1]
    data[2] = data[2].str.split('::').str[0]
    data[4] = data[3].str.split('-').str[1]
    data[3] = data[3].str.split('-').str[0]
    data.columns = ['Sample', 'Chrom', 'Position', 'start', 'end']

    # get dist of mut from start of probe
    data['dist'] = data['Position'].astype(int) - data['start'].astype(int)

    # get probes for each mut where mut_pos is closest to center of probe
    min_dist = data.loc[data.groupby('Position')['dist'].apply(lambda x: (x - size/2).abs().idxmin())]
    
    # reform probe_id names for output file
    selected_ids = min_dist['Sample'] + '_' + min_dist['Chrom'] + '_' + min_dist['Position'] + '::' + min_dist['Chrom'] + ':' + min_dist['start'] + '-' + min_dist['end']
    min_dist['probe_id'] = selected_ids

    # remove old files if they exist and save output files
    path_stats = path_output + "/probe_loc_stats.tsv"
    try:
        os.remove(path_stats)
    except FileNotFoundError:
        pass
    min_dist.to_csv(path_stats, index=False, sep='\t')
    path_min_dist_ids = path_output + "/closest_center_probe_ids.txt"
    try:
        os.remove(path_min_dist_ids)
    except FileNotFoundError:
        pass
    selected_ids.to_csv(path_min_dist_ids, index=False, sep='\t', header=None)

    print(f"List of selected probe_ids: {path_min_dist_ids}")
    print(f"Table of selected probe dist stats: {path_stats}")


def main():
    parser = argparse.ArgumentParser(description="get probes with mut closest to center")
    parser.add_argument("--path_output", required=True, help="abs path to output file")
    parser.add_argument("--path_input", required=True, help="abs path to intput file")
    parser.add_argument("--size", required=False, type=int, default=120, help="length of candidate probes")
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()

    closest_center_probes(args.path_output, args.path_input, args.size)
    
if __name__ == "__main__":
    main()


