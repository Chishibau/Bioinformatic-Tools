#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# File paths
SIGNATURE_FILE = "/groups/wyattgrp/users/zshong/projects/tmp_outputs/bladder_wgs_cluster_order.tsv"
MUTATIONS_FILE = "/groups/wyattgrp/users/zshong/projects/tmp_outputs/mutato_analyze_with_pt.tsv"
OUTPUT_FILE = "/groups/wyattgrp/users/zshong/projects/tmp_outputs/cluster_mut_freq_bar.png"

GENES = ['FGFR3', 'ERBB2', 'PIK3CA', 'PPARG', 'CCND1', 'TACC3', 'TSC1', 'ARID1A', 'KDM6A', 'KMT2D', 'TP53', 'RB1', 'CDKN2A', 'CDKN2B', 'MTAP', 'ERCC2', 'NECTIN4', 'TERT']

# Cluster labels and colors
CLUSTER_INFO = {
    1: {'label': 'Cluster 1', 'color': 'deeppink'},
    2: {'label': 'Cluster 2', 'color': 'black'},
    3: {'label': 'Cluster 3', 'color': 'red'},
    4: {'label': 'Cluster 4', 'color': 'mediumaquamarine'},
    5: {'label': 'Cluster 5', 'color': 'mediumpurple'},
    6: {'label': 'Cluster 6', 'color': 'orange'}
}

def main():
    # Load data
    signatures = pd.read_csv(SIGNATURE_FILE, sep='\t')
    mutations = pd.read_csv(MUTATIONS_FILE, sep='\t')
    
    # Extract patient from sample name
    mutations['Patient'] = mutations['Sample'].str.extract(r'^(GU-\d+-\d+)')[0]
    signatures['Patient'] = signatures['Sample'].str.extract(r'^(GU-\d+-\d+)')[0]
    
    # Check for patients with samples in different clusters
    patient_clusters = signatures.groupby('Patient')['cluster'].nunique()
    patients_to_exclude = patient_clusters[patient_clusters > 1].index.tolist()
    
    print(f"Excluding {len(patients_to_exclude)} patients with samples in different clusters:")
    print(patients_to_exclude)
    
    # Filter out these patients
    mutations = mutations[~mutations['Patient'].isin(patients_to_exclude)]
    signatures = signatures[~signatures['Patient'].isin(patients_to_exclude)]
    
    # Merge with cluster information
    mutations = mutations.merge(signatures[['Sample', 'cluster']], on='Sample', how='left')
    
    # Filter for genes of interest
    mutations = mutations[mutations['Gene'].isin(GENES)]
    
    # Get unique patients per cluster
    n_patients = signatures.groupby('cluster')['Patient'].nunique()
    
    # Count patients with mutations per gene per cluster
    mutation_counts = mutations.groupby(['cluster', 'Gene'])['Patient'].nunique().reset_index()
    mutation_counts.columns = ['cluster', 'Gene', 'n_mutated']
    
    # Calculate frequency as percentage
    mutation_freq = mutation_counts.merge(
        n_patients.reset_index().rename(columns={'Patient': 'n_total'}),
        on='cluster'
    )
    mutation_freq['frequency'] = (mutation_freq['n_mutated'] / mutation_freq['n_total']) * 100
    
    # Create a complete grid of all genes x all clusters
    all_clusters = list(CLUSTER_INFO.keys())
    complete_grid = pd.DataFrame(
        [(g, c) for g in GENES for c in all_clusters],
        columns=['Gene', 'cluster']
    )
    
    # Merge with actual frequencies (missing combinations will be 0)
    freq_complete = complete_grid.merge(
        mutation_freq[['Gene', 'cluster', 'frequency']], 
        on=['Gene', 'cluster'], 
        how='left'
    ).fillna(0)
    
    # Pivot for plotting
    freq_pivot = freq_complete.pivot(index='Gene', columns='cluster', values='frequency')
    
    # Reorder genes to match the gene list
    freq_pivot = freq_pivot.loc[GENES]
    
    # Create bar plot
    fig, ax = plt.subplots(figsize=(16, 7))
    
    # Bar width and positions - add more space between genes
    n_clusters = len(all_clusters)
    x = np.arange(len(GENES)) * 1.5  # Multiply by 1.5 to add space between gene groups
    width = 0.15  # Slightly wider bars since we have more space
    
    # Plot bars for each cluster
    for i, cluster in enumerate(all_clusters):
        offset = width * (i - n_clusters/2 + 0.5)
        color = CLUSTER_INFO[cluster]['color']
        label = CLUSTER_INFO[cluster]['label']
        
        ax.bar(x + offset, freq_pivot[cluster], width, 
               label=label, color=color, edgecolor='white', linewidth=0.5)
    
    # Customize plot
    ax.set_ylabel('Mutation frequency\n(% patients)', fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(GENES, fontsize=11)
    ax.set_ylim(0, 100)
    ax.legend(fontsize=10, frameon=True, loc='upper right', bbox_to_anchor=(1.05, 1))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    plt.title("Bladder WGS mutational frequency by Spearman clusters")
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {OUTPUT_FILE}")
    
    # Print summary
    print(f"\nPatients per cluster:")
    print(n_patients)
    
    plt.show()

if __name__ == "__main__":
    main()