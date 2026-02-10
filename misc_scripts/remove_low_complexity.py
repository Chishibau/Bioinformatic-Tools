#!/usr/bin/env python3

"""
Filter mutations that fall in soft-masked (lowercase) regions of a reference genome.

Lowercase bases in FASTA files typically indicate repetitive or low-complexity regions
that have been soft-masked by RepeatMasker or similar tools.

Usage:
    python /groups/wyattgrp/users/zshong/codebook/misc_scripts/remove_low_complexity.py <mutations_file> <output_file> /groups/wyattgrp/reference/hg38/hg38.fa [options]

python /groups/wyattgrp/users/zshong/codebook/misc_scripts/remove_low_complexity.py /groups/wyattgrp/projects/brachy_tracks/mutations/mut_analyze_bq20/muts_noncoding_VAF.tsv /groups/wyattgrp/projects/brachy_tracks/probe_selection/thresholding/noncoding/noncoding_lowcomp_removed.tsv /groups/wyattgrp/reference/hg38/hg38.fa
    
"""
import pandas as pd
import sys
import argparse


def is_softmasked(chrom, pos, ref_fasta, window=0):
    """
    Check if a genomic position is in a soft-masked (lowercase) region.
    
    Parameters:
    -----------
    chrom : str
        Chromosome name
    pos : int
        Genomic position (1-based)
    ref_fasta : pysam.FastaFile
        Opened reference FASTA file
    window : int
        Check surrounding bases as well (default: 0 = only check exact position)
    
    Returns:
    --------
    bool : True if position (or any base in window) is lowercase
    """
    try:
        # Get the base at the position (pysam uses 0-based coordinates)
        start = max(0, pos - 1 - window)
        end = pos + window
        
        seq = ref_fasta.fetch(chrom, start, end)
        
        # Check if any base in the region is lowercase
        if any(c.islower() for c in seq):
            return True
        
        return False
        
    except Exception as e:
        print(f"Warning: Could not fetch {chrom}:{pos} - {e}")
        return False


def filter_softmasked_mutations(mutations_file, output_file, ref_fasta_path, window=0, keep_flag=False):
    """
    Filter out mutations in soft-masked regions.
    
    Parameters:
    -----------
    mutations_file : str
        Path to input mutations TSV file
    output_file : str
        Path to output filtered TSV file
    ref_fasta_path : str
        Path to reference genome FASTA file
    window : int
        Check surrounding bases (0 = only exact position)
    keep_flag : bool
        Keep the is_softmasked column in output
    """
    
    try:
        import pysam
    except ImportError:
        print("ERROR: pysam is required. Install with: pip install pysam")
        sys.exit(1)
    
    print(f"Reading mutations from: {mutations_file}")
    
    # Read mutations file
    try:
        df = pd.read_csv(mutations_file, sep='\t')
    except:
        df = pd.read_csv(mutations_file, sep=',')
    
    print(f"Loaded {len(df)} mutations")
    
    # Find chromosome and position columns (case-insensitive)
    chrom_col = None
    pos_col = None
    
    for col in df.columns:
        col_lower = col.lower()
        if 'chrom' in col_lower and chrom_col is None:
            chrom_col = col
        elif 'pos' in col_lower and pos_col is None:
            pos_col = col
    
    if chrom_col is None or pos_col is None:
        print(f"ERROR: Could not find chromosome and/or position columns")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)
    
    print(f"Using columns: Chrom={chrom_col}, Pos={pos_col}")
    print(f"Opening reference genome: {ref_fasta_path}")
    
    # Open reference FASTA
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    
    # Check each mutation
    is_softmasked_list = []
    
    for idx, row in df.iterrows():
        if idx % 10000 == 0 and idx > 0:
            print(f"  Processed {idx}/{len(df)} mutations...")
        
        chrom = str(row[chrom_col])
        pos = int(row[pos_col])
        
        # Check if position is soft-masked
        is_masked = is_softmasked(chrom, pos, ref_fasta, window)
        is_softmasked_list.append(is_masked)
    
    ref_fasta.close()
    
    # Add soft-masked flag
    df['is_softmasked'] = is_softmasked_list
    
    # Filter
    filtered_df = df[~df['is_softmasked']].copy()
    
    # Remove flag column unless requested to keep
    if not keep_flag:
        filtered_df = filtered_df.drop('is_softmasked', axis=1)
    
    # Report results
    n_masked = sum(is_softmasked_list)
    
    print(f"\nResults:")
    print(f"  Total mutations: {len(df)}")
    print(f"  In soft-masked regions (lowercase): {n_masked} ({100*n_masked/len(df):.1f}%)")
    print(f"  In non-masked regions (uppercase): {len(filtered_df)} ({100*len(filtered_df)/len(df):.1f}%)")
    print(f"  Removed: {len(df) - len(filtered_df)}")
    
    # Save filtered mutations
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nFiltered mutations saved to: {output_file}")
    
    # Optionally save the masked mutations
    if n_masked > 0:
        masked_file = output_file.replace('.tsv', '_softmasked.tsv')
        df[df['is_softmasked']].drop('is_softmasked', axis=1).to_csv(masked_file, sep='\t', index=False)
        print(f"Soft-masked mutations saved to: {masked_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Filter mutations in soft-masked (lowercase) regions of reference genome'
    )
    
    parser.add_argument('mutations', help='Input mutations file (TSV)')
    parser.add_argument('output', help='Output filtered mutations file (TSV)')
    parser.add_argument('reference', help='Reference genome FASTA file')
    parser.add_argument('--window', type=int, default=0,
                       help='Check surrounding bases (0 = only exact position, default: 0)')
    parser.add_argument('--keep-flag', action='store_true',
                       help='Keep the is_softmasked column in output')
    
    args = parser.parse_args()
    
    filter_softmasked_mutations(
        args.mutations,
        args.output,
        args.reference,
        window=args.window,
        keep_flag=args.keep_flag
    )