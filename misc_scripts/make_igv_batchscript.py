#!/usr/bin/env python3

"""
IGV Batch Script Generator
Creates IGV batch scripts for visualizing mutations from a mutations file.

Usage:
    python generate_igv_batch.py <mutations_file> <output_script> <cram_dir> [options]
    
Example:
    python generate_igv_batch.py mutations.tsv igv_batch.txt /path/to/crams
    python generate_igv_batch.py mutations.tsv igv_batch.txt /path/to/crams --snapshot-dir /path/to/snapshots --wbc-map samples.tsv

python /groups/wyattgrp/users/zshong/codebook/misc_scripts/make_igv_batchscript.py /groups/wyattgrp/projects/brachy_tracks/mutations/igv_coding_VAF80p_depth20p/rerun_muts.tsv /groups/wyattgrp/projects/brachy_tracks/mutations/igv_coding_VAF80p_depth20p/igv_batchscript.txt /groups/wyattgrp/projects/brachy_tracks/cram_alignments --snapshot-dir /groups/wyattgrp/projects/brachy_tracks/mutations/igv_coding_VAF80p_depth20p/ --wbc-map /groups/wyattgrp/projects/brachy_tracks/reference/sample_wbc_mapping.tsv --genome /groups/wyattgrp/reference/hg38/hg38.fa
"""

import pandas as pd
import sys
import os
from pathlib import Path


def generate_igv_batch_script(mutations_file, output_script, cram_dir, config=None):
    """
    Generate an IGV batch script from a mutations file.
    
    Parameters:
    -----------
    mutations_file : str
        Path to TSV file containing mutations with columns:
        Sample, Chrom, Position, Ref allele, Gene, Effect
    output_script : str
        Path where the IGV batch script will be written
    cram_dir : str
        Directory containing CRAM files (REQUIRED)
    config : dict, optional
        Configuration dictionary with keys:
        - snapshot_dir: Directory for IGV snapshots (default: same as output_script directory)
        - genome: Reference genome (default: hg38)
        - sample_wbc_map: Dictionary mapping sample names to WBC names (default: None)
    """
    
    # Default configuration
    default_config = {
        'snapshot_dir': None,  # Will default to output script directory
        'genome': 'hg38',
        'sample_wbc_map': None,
    }
    
    if config is None:
        config = {}
    
    # Merge with defaults
    config = {**default_config, **config}
    
    # Validate cram_dir exists
    if not os.path.isdir(cram_dir):
        print(f"ERROR: CRAM directory not found: {cram_dir}")
        sys.exit(1)
    
    # Set snapshot directory to output directory if not specified
    if config['snapshot_dir'] is None:
        config['snapshot_dir'] = str(Path(output_script).parent.absolute())
    
    print(f"Reading mutations from: {mutations_file}")
    print(f"CRAM directory: {cram_dir}")
    
    # Read mutations file - try different separators
    try:
        muts = pd.read_csv(mutations_file, sep='\t')
    except:
        try:
            muts = pd.read_csv(mutations_file, sep=',')
        except Exception as e:
            print(f"Error reading file: {e}")
            raise
    
    print(f"Loaded {len(muts)} mutations")
    
    # Validate required columns
    required_cols = ['Sample', 'Chrom', 'Position', 'Ref allele']
    missing_cols = [col for col in required_cols if col not in muts.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        print(f"Available columns: {list(muts.columns)}")
        sys.exit(1)
    
    # Optional columns with defaults
    if 'Gene' not in muts.columns:
        print("Warning: 'Gene' column not found, using 'Unknown' for all mutations")
        muts['Gene'] = 'Unknown'
    
    if 'Effect' not in muts.columns:
        print("Warning: 'Effect' column not found, using 'Unknown' for all mutations")
        muts['Effect'] = 'Unknown'
    
    # Clean up data
    muts['Gene'] = muts['Gene'].fillna('Intergenic')
    muts['Effect'] = muts['Effect'].fillna('Unknown')
    
    # Filter to valid chromosomes if Chrom contains chromosome info
    if muts['Chrom'].dtype == 'object':
        original_len = len(muts)
        muts = muts[muts['Chrom'].astype(str).str.contains('chr', case=False, na=False)]
        filtered = original_len - len(muts)
        if filtered > 0:
            print(f"Filtered out {filtered} mutations on non-standard chromosomes")
    
    # Filter out rows without probe_id if that column exists
    if 'probe_id' in muts.columns:
        original_len = len(muts)
        muts = muts[~muts['probe_id'].isna()]
        filtered = original_len - len(muts)
        if filtered > 0:
            print(f"Filtered out {filtered} mutations without probe_id")
    
    print(f"Generating IGV batch script for {len(muts)} mutations across {muts['Sample'].nunique()} samples")
    
    # Create output directory if it doesn't exist
    os.makedirs(config['snapshot_dir'], exist_ok=True)
    
    # Write IGV batch script
    with open(output_script, 'w') as f:
        
        for sample, sample_muts in muts.groupby('Sample'):
            print(f"  Processing sample: {sample} ({len(sample_muts)} mutations)")
            
            # Path to sample CRAM file
            cram_path = os.path.join(cram_dir, f"{sample}.cram")
            crai_path = f"{cram_path}.crai"
            
            # Check if CRAM file exists
            if not os.path.exists(cram_path):
                print(f"    WARNING: CRAM file not found: {cram_path}")
            
            # Start new session
            f.write('new\n')
            
            # Load sample CRAM
            f.write(f"load {cram_path} index={crai_path}\n")
            
            # Load WBC CRAM if mapping is provided
            if config['sample_wbc_map'] and sample in config['sample_wbc_map']:
                wbc_sample = config['sample_wbc_map'][sample]
                wbc_cram = os.path.join(cram_dir, f"{wbc_sample}.cram")
                wbc_crai = f"{wbc_cram}.crai"
                f.write(f"load {wbc_cram} index={wbc_crai}\n")
            
            # Set snapshot directory and genome
            f.write(f'snapshotDirectory {config["snapshot_dir"]}\n')
            f.write(f'genome {config["genome"]}\n')
            
            # Process each mutation
            for idx, mut in sample_muts.iterrows():
                # Extract chromosome number (remove 'chr' prefix if present)
                chrom = str(mut['Chrom'])
                if chrom.lower().startswith('chr'):
                    chrom_num = chrom[3:]
                else:
                    chrom_num = chrom
                
                # Calculate IGV position (adjust for ref allele length)
                real_pos = int(mut['Position'])
                ref_allele = str(mut['Ref allele'])
                igv_pos = real_pos + len(ref_allele) - 1
                
                # Go to position
                f.write(f'goto {chrom_num}:{igv_pos}\n')
                f.write('sort base\n')
                f.write('colorBy READ_STRAND\n')
                
                # Create snapshot filename
                gene = str(mut['Gene']).replace(' ', '_')
                effect = str(mut['Effect']).split(' ')[0]
                snapshot_name = f"{sample},{gene}.{chrom_num}.{real_pos}.{effect}.png"
                
                f.write(f'snapshot {snapshot_name}\n')

        # Close IGV after all screenshots are saved
        f.write('\nexit\n')
        
    print(f"\nIGV batch script written to: {output_script}")
    print(f"Snapshots will be saved to: {config['snapshot_dir']}")


def load_sample_wbc_mapping(mapping_file):
    """
    Load sample to WBC mapping from a file.
    
    Parameters:
    -----------
    mapping_file : str
        Path to TSV/CSV file with 'sample' and 'wbc' columns
        
    Returns:
    --------
    dict : Dictionary mapping sample names to WBC names
    """
    try:
        df = pd.read_csv(mapping_file, sep='\t')
    except:
        df = pd.read_csv(mapping_file, sep=',')
    
    if 'sample' not in df.columns or 'wbc' not in df.columns:
        print(f"Warning: mapping file should have 'sample' and 'wbc' columns")
        return None
    
    return dict(zip(df['sample'], df['wbc']))


if __name__ == "__main__":
    
    if len(sys.argv) < 4:
        print("Usage: python generate_igv_batch.py <mutations_file> <output_script> <cram_dir> [options]")
        print("\nRequired Arguments:")
        print("  mutations_file           TSV file with mutations")
        print("  output_script            Path for output IGV batch script")
        print("  cram_dir                 Directory containing CRAM files")
        print("\nOptional Arguments:")
        print("  --snapshot-dir <path>    Directory for snapshots (default: same as output)")
        print("  --genome <name>          Reference genome (default: hg38)")
        print("  --wbc-map <file>         TSV file with sample-to-WBC mapping")
        print("\nExamples:")
        print("  python generate_igv_batch.py mutations.tsv igv_batch.txt /path/to/crams")
        print("  python generate_igv_batch.py mutations.tsv igv_batch.txt /path/to/crams --snapshot-dir /path/to/snapshots")
        print("  python generate_igv_batch.py mutations.tsv igv_batch.txt /path/to/crams --wbc-map samples.tsv --genome hg19")
        sys.exit(1)
    
    mutations_file = sys.argv[1]
    output_script = sys.argv[2]
    cram_dir = sys.argv[3]
    
    # Parse optional arguments
    config = {}
    i = 4
    while i < len(sys.argv):
        if sys.argv[i] == '--snapshot-dir' and i + 1 < len(sys.argv):
            config['snapshot_dir'] = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--genome' and i + 1 < len(sys.argv):
            config['genome'] = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--wbc-map' and i + 1 < len(sys.argv):
            config['sample_wbc_map'] = load_sample_wbc_mapping(sys.argv[i + 1])
            i += 2
        else:
            print(f"Unknown argument: {sys.argv[i]}")
            i += 1
    
    # Validate input file exists
    if not os.path.exists(mutations_file):
        print(f"ERROR: Mutations file not found: {mutations_file}")
        sys.exit(1)
    
    # Generate the script
    generate_igv_batch_script(mutations_file, output_script, cram_dir, config)