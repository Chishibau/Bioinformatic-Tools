# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 14:00:38 2025

@author: zshong
"""

""" Example usage
python /groups/wyattgrp/users/zshong/codebook/plotting_scripts/fix_seg_error_allgenes.py \
    --fix_seg_dir /groups/wyattgrp/users/amunzur/hla_project/copynum/fix/genes \
    --out_path /groups/wyattgrp/users/zshong/codebook/copynum_violin

NOTE: automatically ignores HLA-A, HLA-B, and HLA-C will need to use separate script
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse
import math

def plot(root, out_path):
    # find all fix files within directory
    files = os.listdir(root)
    fix_files = list(filter(lambda file: file.endswith(".fix"), files))

    # %%

    # initialize sample value dfs
    fix = pd.DataFrame()
    seg = pd.DataFrame()

    # calculate variation from segments for each sample
    error_dict = {} # init dict to store error values
    sample_count_dict = {} # init dict for num of samples per pos per gene
    exclude_genes = ("HLA-A", "HLA-B", "HLA-C") # exclude HLA genes

    for sample in fix_files:
        fix = pd.read_csv(os.path.join(root, sample), sep='\t')
        seg = pd.read_csv(os.path.join(root, sample+".segment"), sep='\t')
        for row in fix.itertuples(index=False):
            # treat intronic muts the same as exonic
            gene = row.gene
            if 'intronic' in gene.lower():
                gene = row.gene.split(' ')[0]
            # ignore hla_a, hla_b, and hla_c muts
            if gene in exclude_genes:
                continue
            # get median position for each bin in fix file
            pos = (row.start + row.end) / 2
            for seg_row in seg.itertuples(index=False):
                if row.start < seg_row.end:
                    # calculate difference between log2 in fix file and corresponding segment
                    err = abs(row.log2 - seg_row.log2)
                    break
                
            # update dictionaries
            if gene not in error_dict:
                error_dict[gene] = {}
                sample_count_dict[gene] = {}

            error_dict[gene][pos] = error_dict[gene].get(pos, 0) + err
            sample_count_dict[gene][pos] = sample_count_dict[gene].get(pos, 0) + 1
        
    # setup subplots for each allele
    n_genes = len(error_dict)
    cols = 5    # 5 plots per row
    rows = math.ceil(n_genes / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows), squeeze=False, sharey=True)
    # flatten axis array for looping
    axes = axes.ravel()

    # calculate average error for each gene across all positions and plot
    for idx, gene in enumerate(error_dict):
        error_dict[gene] = {
        pos: (error_dict[gene][pos] / sample_count_dict[gene][pos])
        for pos in error_dict[gene]
        }
        ax = axes[idx]
        ax.set_title(f'{gene}, positions={len(error_dict[gene])}', fontsize=10)    
        ax.plot(list(error_dict[gene].keys()), list(error_dict[gene].values()), '-', linewidth=1)

    for ax in axes[n_genes:]:
        ax.axis("off")

    total_n = len(fix_files)
    fig.suptitle(f'All Genes, n={total_n}\nAverage absolute difference between cnvkit fix log2 values and corresponding segment\nx=bp location within genome, y=avg log2 value differenece\n{root}', 
                fontsize=12, y=1.01)

    plt.subplots_adjust(hspace=0.2)
    plt.tight_layout()
    
    fig.savefig(os.path.join(out_path, "all_genes_fix_error") + '.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(out_path, "all_genes_fix_error") + '.png', dpi=300, bbox_inches='tight')
    print(f'Plots saved to: {os.path.join(out_path, "all_genes_fix_error") + '.pdf'} and {os.path.join(out_path, "all_genes_fix_error") + '.png'}')

def main():
    parser = argparse.ArgumentParser(
        description="Generate error plots from cnvkit fix and segment files for all genes across all samples within given directory"
    )

    parser.add_argument("--fix_seg_dir", required=True, help="absolute path to dir containing cnvkit fix and segmentation files")
    parser.add_argument("--out_path", required=True, help="absolute path to directory for plot")

    args = parser.parse_args()
    
    os.makedirs(args.out_path, exist_ok=True)

    plot(args.fix_seg_dir, args.out_path)

if __name__ == "__main__":
    main()