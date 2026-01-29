#!/usr/bin/env python3

"""
Created on Dec 6 2023
@author: mstephenson, amurtha
modified: zshong Jan 2026
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import sys

# =============================================================================
# CONFIGURATION - Edit these variables
# =============================================================================
SAMPLE_ORDER_FILE = "/groups/wyattgrp/users/zshong/projects/tmp_outputs/wgs_bladder_sample_order.txt"  # Set to path like "/path/to/sample_order.txt" or None to auto-sort
INPUT_FILE = "/groups/wyattgrp/users/zshong/projects/tmp_outputs/mutato_analyze_with_pt.tsv"
OUTPUT_FILE = "/groups/wyattgrp/users/zshong/projects/tmp_outputs/oncogone_bladder.png"
CLUSTER_ANNOTATION_FILE = "/groups/wyattgrp/users/zshong/projects/tmp_outputs/bladder_wgs_cluster_order.tsv"  # Set to path like "/path/to/clusters.tsv" or None to skip

def keep_coding_muts(df):
    coding = ['Stopgain',  'Missense', 'Splice',
           'Frameshift',  'Non-frameshift', "Promoter"]
    non_coding = [ "3'-UTR", 'Intronic',  'Upstream','Downstream', 'Synonymous', "5'-UTR"]
    print("Missing effects", [e for e in df['Effect'].str.split(' ').str[0].unique().tolist() if e not in coding+non_coding])
    df = df[df['Effect'].str.split(' ').str[0].isin(coding)]
    return df

def main():
    muts_color = {'Missense':'#79B443',
              'Frameshift':'#FFC907',
              'Stopgain':'#FFC907',
              'Non-frameshift':'#BD4398',
              'Splice':'#FFC907',
              'Promoter':'#FF3300'}
    
    genes = ['FGFR3', 'ERBB2', 'PIK3CA', 'PPARG', 'CCND1', 'TACC3', 'TSC1', 'ARID1A', 'KDM6A', 'KMT2D', 'TP53', 'RB1', 'CDKN2A', 'CDKN2B', 'MTAP', 'ERCC2', 'NECTIN4', 'TERT']
    g_order = dict(zip(genes, np.arange(len(genes))))
    
    muts = pd.read_csv(INPUT_FILE, sep='\t')
    muts = muts[muts["Sample"].str.contains("cfDNA", na=False)]
    
    # Keep coding mutations in genes first
    muts = muts[muts['Gene'].isin(genes)]
    muts = keep_coding_muts(muts)
    
    # Determine sample order
    if SAMPLE_ORDER_FILE:
        # Read sample order from file (one sample per line)
        with open(SAMPLE_ORDER_FILE, 'r') as f:
            sample_order_list = [line.strip() for line in f if line.strip()]
        
        # Filter to only samples that exist in data
        existing_samples = set(muts['Sample'].unique())
        sample_order_list = [s for s in sample_order_list if s in existing_samples]
        
        # Add any samples from data not in the order file
        missing_samples = sorted(list(existing_samples - set(sample_order_list)))
        if missing_samples:
            print(f"Warning: {len(missing_samples)} samples in data but not in order file. Adding to end.")
            sample_order_list.extend(missing_samples)
    else:
        # Auto-sort: Order samples by patient, then by sample name
        muts_sorted = muts.sort_values(['Patient', 'Sample'])
        sample_order_list = muts_sorted['Sample'].unique().tolist()
    
    s_order = dict(zip(sample_order_list, range(0, len(sample_order_list))))
    
    muts['x'] = muts['Sample'].map(s_order)
    muts['y'] = muts['Gene'].map(g_order)
    
    # Adjust X and Y for multiple mutations - deduplicate by color within sample/gene
    # First, keep only unique color per sample/gene combination
    muts['color'] = muts['Effect'].str.split(' ').str[0].map(muts_color)
    muts_dedup = muts.groupby(['Sample', 'Gene', 'color']).first().reset_index()
    
    # Now count unique colors per sample/gene
    m_count_pt = muts_dedup.groupby(['Sample','Gene']).size().reset_index(name='gene_mut_count')
    muts_dedup = muts_dedup.merge(m_count_pt, on=['Sample','Gene'], how='left')
    
    offset = 0.2
    
    for index, row in m_count_pt.iterrows():
        if row['gene_mut_count'] > 1:
            gene = row['Gene']
            s = row['Sample']
            colors = muts_dedup[(muts_dedup['Gene'] == gene)&(muts_dedup['Sample'] == s)]['color'].to_list()
            for i, color in enumerate(colors):
                if i == 0:
                    muts_dedup.loc[(muts_dedup['Gene'] == gene) & (muts_dedup['color'] == color) & (muts_dedup['Sample'] == s), 'y'] = muts_dedup['y']-offset
                elif i == 1:
                    muts_dedup.loc[(muts_dedup['Gene'] == gene) & (muts_dedup['color'] == color) & (muts_dedup['Sample'] == s), 'y'] = muts_dedup['y']+offset
                else:
                    # For 3+ mutations, stack vertically
                    muts_dedup.loc[(muts_dedup['Gene'] == gene) & (muts_dedup['color'] == color) & (muts_dedup['Sample'] == s), 'y'] = muts_dedup['y'] + (i-1)*0.15
    
    # Use deduplicated data for plotting
    muts = muts_dedup
    
    # Load cluster annotations if provided
    cluster_assignment = None
    if CLUSTER_ANNOTATION_FILE:
        clusters = pd.read_csv(CLUSTER_ANNOTATION_FILE, sep='\t')
        # Map cluster colors and assignments to samples
        cluster_map = dict(zip(clusters['Sample'], clusters['cluster_color']))
        cluster_assignment = dict(zip(clusters['Sample'], clusters['cluster']))
        
        # If no custom sample order file, use order from cluster file
        if not SAMPLE_ORDER_FILE:
            cluster_sample_order = clusters['Sample'].tolist()
            # Filter to only samples that exist in mutation data
            existing_samples = set(muts['Sample'].unique())
            sample_order_list = [s for s in cluster_sample_order if s in existing_samples]
            
            # Add any samples from data not in cluster file
            missing_samples = sorted(list(existing_samples - set(sample_order_list)))
            if missing_samples:
                print(f"Warning: {len(missing_samples)} samples in mutation data but not in cluster file. Adding to end.")
                sample_order_list.extend(missing_samples)
            
            s_order = dict(zip(sample_order_list, range(0, len(sample_order_list))))
            muts['x'] = muts['Sample'].map(s_order)
    else:
        clusters = None
        cluster_map = None
    
    # Create figure with annotation track if clusters provided
    if clusters is not None:
        fig, (ax_top, ax_bottom) = plt.subplots(nrows=2, figsize=(14, 14), 
                                                 gridspec_kw={'height_ratios':[0.3, 6]},
                                                 sharex=True)
    else:
        fig, ax_bottom = plt.subplots(figsize=(14, 14))
        ax_top = None
    
    # Add grey background boxes for each gene x sample position
    for i in range(len(s_order)):
        for j in range(len(genes)):
            rect = plt.Rectangle((i-0.4, j-0.4), 0.8, 0.8, 
                                facecolor='#E6E7E8', edgecolor='white', 
                                linewidth=1, zorder=1)
            ax_bottom.add_patch(rect)
    
    # Plot mutations on top - smaller markers
    ax_bottom.scatter(muts['x'], muts['y'], c=muts['color'], marker='s', s=40, 
              edgecolors='white', linewidths=0.5, zorder=10)
    
    # Set up x-axis with sample names - align properly
    sample_names = [s for s in s_order.keys()]
    ax_bottom.set_xticks(range(len(sample_names)))
    ax_bottom.set_xticklabels(sample_names, rotation=90, ha='center', va='top', fontsize=7)
    # ax_bottom.tick_params(length='')  # Remove tick marks but keep labels
    ax_bottom.set_xlim(-0.5, len(sample_names)-0.5)  # Ensure proper alignment
    
    # Set up y-axis with gene names
    ax_bottom.set_ylim(-0.5, len(genes)-0.5)
    ax_bottom.set_yticks(list(g_order.values()))
    ax_bottom.set_yticklabels(list(g_order.keys()), fontsize=11)
    ax_bottom.tick_params(left=False, which='both')  # Remove tick marks but keep labels
    ax_bottom.invert_yaxis()
    
    # Add gridlines
    ax_bottom.set_axisbelow(True)
    
    # Clean up spines
    ax_bottom.spines['top'].set_visible(False)
    ax_bottom.spines['right'].set_visible(False)
    ax_bottom.spines['left'].set_visible(False)
    ax_bottom.spines['bottom'].set_visible(False)
    
    # Add vertical lines between clusters
    if cluster_assignment is not None:
        prev_cluster = None
        for i, sample in enumerate(sample_names):
            curr_cluster = cluster_assignment.get(sample)
            if prev_cluster is not None and curr_cluster != prev_cluster:
                # Draw vertical line between clusters
                ax_bottom.axvline(x=i-0.5, color='black', linewidth=1, linestyle='--', zorder=15)
            prev_cluster = curr_cluster
    
    # Plot cluster annotation track if provided
    if ax_top is not None and cluster_map is not None:
        for i, sample in enumerate(sample_names):
            color = cluster_map.get(sample, '#CCCCCC')  # Default gray if not found
            rect = plt.Rectangle((i-0.4, -0.4), 0.8, 0.8, 
                                facecolor=color, edgecolor='white', 
                                linewidth=1)
            ax_top.add_patch(rect)
        
        ax_top.set_ylim(-0.5, 0.5)
        ax_top.set_xlim(-1, len(sample_names)-0.5)
        ax_top.set_yticks([0])
        ax_top.set_yticklabels(['Group'], fontsize=10)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        ax_top.spines['left'].set_visible(False)
        ax_top.spines['bottom'].set_visible(False)
        ax_top.tick_params(left=False, bottom=False, top=False)
    
    # Add legend - position outside plot area
    def l(c, label, m='s', ms=8):
        return Line2D([], [], lw=0, color=c, marker=m, markeredgewidth=0.5, 
                     markeredgecolor='white', ms=ms, label=label)
    
    legend_elements = [
        l('#79B443', 'Missense'),
        l('#FFC907', 'Truncating (Stopgain/Frameshift/Splice)'),
        l('#BD4398', 'In-frame InDel'),
        l('#FF3300', 'Promoter')
    ]
    ax_bottom.legend(handles=legend_elements, bbox_to_anchor=(1, 0.5), 
                    frameon=True, fontsize=7)
    
    ax_bottom.set_xlabel('Samples', fontsize=11)
    ax_bottom.set_ylabel('Genes', fontsize=11)
    
    # Add title at the very top of the figure
    fig.suptitle('Bladder WGS oncoplot based on Spearman cluster groups', fontsize=13, y=0.98)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=200, bbox_inches='tight')
    print(f"Plot saved to: {OUTPUT_FILE}")
    plt.show()

if __name__ == "__main__":
    main()