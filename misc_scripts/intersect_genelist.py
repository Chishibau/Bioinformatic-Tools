#!/usr/bin/env python3

import pandas as pd

file = "/groups/wyattgrp/projects/brachy_tracks/mutations/mut_analyze_mapq20/muts_coding_VAF.tsv"
genes = "/groups/wyattgrp/projects/brachy_tracks/reference/crpc2022_genes_CHIP_excluded.txt"
output = "/groups/wyattgrp/projects/brachy_tracks/probe_selection/target_genes/mapq20/crpc2022_coding_mapq20.tsv"

df = pd.read_csv(file, sep="\t")
with open(genes, 'r') as f:
    gene_list = [line.rstrip() for line in f]
df["read_depth"] = pd.to_numeric(df["read_depth"], errors="coerce")
df["VAF"] = pd.to_numeric(df["VAF"], errors="coerce")
print(gene_list)

df = df[df['Gene'].isin(gene_list)]

df.to_csv(output, sep='\t', index=False)