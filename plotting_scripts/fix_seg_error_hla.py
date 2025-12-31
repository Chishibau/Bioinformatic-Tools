# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 10:28:47 2025

@author: zshong
"""

""" Example usage
python /groups/wyattgrp/users/zshong/codebook/plotting_scripts/fix_seg_error_hla.py \
    --fix_seg_dir /groups/wyattgrp/users/amunzur/hla_project/allele_specific_copynum/fix \
    --hla_gene hla_b \
    --out_path /groups/wyattgrp/users/zshong/codebook/copynum_violin
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import re
import math
import argparse

def plot(root, gene_name, out_path):
    # find all fix/segment files for target gene
    files = os.listdir(root)
    target_gene_files = list(filter(lambda file: gene_name in file, files))
    fix_files = list(filter(lambda file: file.endswith(".fix"), target_gene_files))

    # for HLA gene, find all unique allele variants
    alleles = []
    pattern = fr"{gene_name}_\d{{2}}_\d{{2}}"
    for file in target_gene_files:
        allele = re.search(pattern, file).group()
        alleles.append(allele)

    alleles = list(set(alleles))

    # setup subplots for each allele
    n_alleles = len(alleles)
    cols = 3    # 3 plots per row
    rows = math.ceil(n_alleles / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4 * rows), squeeze=False, sharey=True)

    # flatten axis array for looping
    axes = axes.ravel()

    # initialize sample value dfs
    fix = pd.DataFrame()
    seg = pd.DataFrame()

    # calculate variation from segments for each allele within each sample
    for idx, allele in enumerate(alleles):
        error_data = {} # initialize dictionary for error values
        sample_n = 0 # initialize counter for sample number
        # get all fix files for a specific allele
        allele_fix_files = list(filter(lambda file: allele in file, fix_files))
        for sample in allele_fix_files:
            sample_n += 1
            fix = pd.read_csv(os.path.join(root, sample), sep='\t')
            seg = pd.read_csv(os.path.join(root, sample+".segment"), sep='\t')
            for row in fix.itertuples(index=False):
                # get median position for each bin in fix file
                pos = (row.start + row.end) / 2
                for seg_row in seg.itertuples(index=False):
                    if row.start < seg_row.end:
                        # calculate difference between log2 in fix file and corresponding segment
                        err = abs(row.log2 - seg_row.log2)
                        break
                error_data[pos] = error_data.get(pos, 0) + err
        
        # calculate average error for each position across samples for all alleles
        error_avg = {key: value / sample_n for key, value in error_data.items()}

        # generate line plot
        ax = axes[idx]
        ax.set_title(f'{allele}, n={sample_n}', fontsize=12)    
        ax.plot(list(error_avg.keys()), list(error_avg.values()), '-', linewidth=1)

    # hide unused axis
    for ax in axes[n_alleles:]:
        ax.axis("off")

    total_n = len(fix_files)
    fig.suptitle(f'{gene_name}, n={total_n}\nAverage absolute difference between cnvkit fix log2 values and corresponding segment\nx=bp location within gene, y=avg log2 value differenece\n{root}', 
                fontsize=15)
    plt.tight_layout()

    fig.savefig(os.path.join(out_path, f"{gene_name}_fix_error") + '.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(out_path, f"{gene_name}_fix_error") + '.png', dpi=300, bbox_inches='tight')
    print(f'Plots saved to: {os.path.join(out_path, f"{gene_name}_fix_error") + '.pdf'} and {os.path.join(out_path, f"{gene_name}_fix_error") + '.png'}')

def main():
    parser = argparse.ArgumentParser(
        description="Generate error plots from cnvkit fix and segment files for HLA genes with different alleles"
    )

    parser.add_argument("--fix_seg_dir", required=True, help="absolute path to dir containing cnvkit fix and segmentation files")
    parser.add_argument("--hla_gene", required=True, help="name for target hla gene")
    parser.add_argument("--out_path", required=True, help="absolute path to directory for plots")

    args = parser.parse_args()
    
    os.makedirs(args.out_path, exist_ok=True)

    plot(args.fix_seg_dir, args.hla_gene, args.out_path)


if __name__ == "__main__":
    main()
