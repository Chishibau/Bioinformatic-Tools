#!/usr/bin/env python3
"""
Memory-efficient version for filtering large files.
Output lines from file1 that are not in file2 based on the first 7 columns.


python /groups/wyattgrp/users/zshong/codebook/misc_scripts/output_unique_lines.py file1.tsv file2.tsv output.tsv
"""

import pandas as pd
import sys

def filter_unique_lines_efficient(file1, file2, output_file, sep='\t', chunksize=100000):
    """
    Memory-efficient version using chunked reading for large files.
    
    Parameters:
    -----------
    file1 : str
        Path to first input file
    file2 : str
        Path to second input file (reference for filtering)
    output_file : str
        Path to output file
    sep : str
        Delimiter for files (default: tab)
    chunksize : int
        Number of rows to read at a time (default: 100000)
    """
    
    print(f"Reading {file2} to build reference set...")
    df2 = pd.read_csv(file2, sep=sep)
    
    # Get first 8 column names
    cols_to_compare = df2.columns[:7].tolist()
    print(f"Comparing based on columns: {cols_to_compare}")
    
    # Create a set of tuples from file2's first 8 columns
    file2_keys = set(df2[cols_to_compare].apply(tuple, axis=1))
    print(f"Number of unique keys in file2: {len(file2_keys)}")
    
    # Free memory
    del df2
    
    # Process file1 in chunks
    print(f"\nProcessing {file1} in chunks...")
    first_chunk = True
    total_lines = 0
    unique_lines = 0
    
    for chunk_num, chunk in enumerate(pd.read_csv(file1, sep=sep, chunksize=chunksize), 1):
        total_lines += len(chunk)
        
        # Get column names from first chunk
        if first_chunk:
            cols_to_compare = chunk.columns[:8].tolist()
            first_chunk = False
        
        # Filter chunk
        mask = ~chunk[cols_to_compare].apply(tuple, axis=1).isin(file2_keys)
        unique_chunk = chunk[mask]
        unique_lines += len(unique_chunk)
        
        # Write to file (append mode after first chunk)
        mode = 'w' if chunk_num == 1 else 'a'
        header = chunk_num == 1
        unique_chunk.to_csv(output_file, sep=sep, index=False, mode=mode, header=header)
        
        print(f"  Chunk {chunk_num}: {len(chunk)} lines, {len(unique_chunk)} unique")
    
    print(f"\nSummary:")
    print(f"  Total lines in file1: {total_lines}")
    print(f"  Lines unique to file1: {unique_lines}")
    print(f"  Saved to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) == 4:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
        output = sys.argv[3]
        filter_unique_lines_efficient(file1, file2, output)
    else:
        print("Usage: python /groups/wyattgrp/users/zshong/codebook/misc_scripts/output_unique_lines.py file1.tsv file2.tsv output.tsv")
        # print("\nThis version is optimized for large files using chunked reading.")