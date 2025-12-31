# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 10:28:47 2025

@author: zshong
"""

""" Example usage
python /groups/wyattgrp/users/zshong/codebook/plotting_scripts/signature_plot.py \

"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import re
import math
import argparse

sig_dict = {
    'Aging': ['SBS1', 'SBS5'],

    'MMR': [
        'SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44',
        'ID7',
        'DBS7', 'DBS10'
    ],

    'POL': [
        'SBS10a', 'SBS10b', 'SBS10c', 'SBS10d', 'SBS28',
        'DBS3'
    ],

    'HR': [
        'SBS3',
        'ID6',
        'DBS13'
    ],

    'BER': ['SBS30', 'SBS36'],

    'Chemotherapy': [
        'SBS11', 'SBS25', 'SBS31', 'SBS35',
        'SBS86', 'SBS87', 'SBS90', 'SBS99',
        'DBS5'
    ],

    'Immunosuppressant': ['SBS32'],

    'Treatment': [
        'SBS11', 'SBS25', 'SBS31', 'SBS32', 'SBS35',
        'SBS86', 'SBS87', 'SBS90', 'SBS99'
    ],

    'APOBEC': ['SBS2', 'SBS13'],

    'Tobacco': [
        'SBS4', 'SBS29', 'SBS92',
        'ID3',
        'DBS2'
    ],

    'UV': [
        'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS38',
        'ID13',
        'DBS1'
    ],

    'AA': [
        'SBS22a', 'SBS22b',
        'ID23',
        'DBS20'
    ],

    'Colibactin': [
        'SBS88',
        'ID18'
    ],

    'Artifact': [
        'SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49',
        'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56',
        'SBS57', 'SBS58', 'SBS59', 'SBS60', 'SBS95',
        'DBS14'
    ],

    'Other': [
        'SBS8', 'SBS12', 'SBS16', 'SBS17a', 'SBS17b',
        'SBS18', 'SBS19', 'SBS23', 'SBS24', 'SBS33', 'SBS34', 'SBS37',
        'SBS39', 'SBS40a', 'SBS40b', 'SBS40c', 'SBS41', 'SBS42',
        'SBS89', 'SBS91', 'SBS93', 'SBS94', 'SBS96', 'SBS97', 'SBS98',
        'ID4', 'ID5', 'ID8', 'ID9', 'ID10', 'ID11', 'ID12', 'ID14',
        'ID15', 'ID16', 'ID19', 'ID20', 'ID21', 'ID22',
        'DBS4', 'DBS6', 'DBS8', 'DBS9', 'DBS11', 'DBS12',
        'DBS15', 'DBS16', 'DBS17', 'DBS18', 'DBS19',
        'DBS21', 'DBS22'
    ],

    'Slippage': ['ID1', 'ID2'],

    'TOP2A': ['ID17']
}


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
