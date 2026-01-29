#!/usr/bin/env python3

"""
Created on Fri Nov 14 11:12:46 2025

@author: zshong
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
import argparse

def parse_dbNSFPS(out_dir, tsv_dir):
    line_counts = []
    
    for file_path in Path(tsv_dir).glob("*.tsv"):
        df = pd.read_csv(file_path, sep="\t")
        original_count = len(df)
        
        avg = df.iloc[:, np.r_[58,90,118,135,142,164]]
        avg = avg.replace('.', np.nan).astype(float)
        
        # Count non-NaN values per row
        non_nan_count = avg.notna().sum(axis=1)
        
        # Filter: keep only rows with at least 2 numeric values
        mask = non_nan_count >= 2
        df = df[mask]
        avg = avg[mask]
        
        filtered_count = len(df)
        
        # Calculate average
        avg = avg.mean(axis=1, numeric_only=True).round(5)
        
        df = df.iloc[:, np.r_[0:7,58,90,118,135,142,164]]
        df.iloc[:,0] = "chr" + df.iloc[:,0].astype(str)
        df['avg'] = avg.values
        df.insert(loc=0, column='Sample', value=file_path.stem)
        df.to_csv(os.path.join(out_dir, f"{file_path.stem}_dbNSFPS.tsv"), index=False, sep="\t")
        
        # Record counts
        line_counts.append({
            'Sample': file_path.stem,
            'Original': original_count,
            'Filtered': filtered_count,
            'Removed': original_count - filtered_count
        })
    
    # Write line counts to file
    counts_df = pd.DataFrame(line_counts)
    counts_df.to_csv(os.path.join(out_dir, "summary_counts.txt"), sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(
        description="parses dbNSFPS output files and thresholds to specified average value for mutation effect scores"
    )

    parser.add_argument("--tsv_dir", required=True, help="Directory containing dbNSFPS output files")
    parser.add_argument("--out_dir", required=True, help="Output directory for parsed dbNSFPS output files")

    args = parser.parse_args()

    # check if the given output directory exists 
    if os.path.isdir(args.out_dir):
        parse_dbNSFPS(args.out_dir, args.tsv_dir)
    else:
        args.out_dir.mkdir(parents=True, exist_ok=True)

if __name__ == "__main__":
    main()