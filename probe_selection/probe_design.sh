#!/bin/bash

# generate 200 bp regions around each mutation postion
python /groups/wyattgrp/users/zshong/codebook/probe_design/bed_from_mut_pos.py \
    --path_output /groups/wyattgrp/projects/chordoma/panel_design/candidates \
    --path_tsv /groups/wyattgrp/projects/chordoma/panel_design/candidates/all_muts.tsv \
    --chrom_col "CHROM" \ # optional, column name of chrosomosome label
    --pos_col "POSITION" \ # optional, column name of mut pos
    --id_col "patient" \ # optional, column name of sample identifier
    --size 100 \ # optional, num of bp to include before/after mut pos

# make 120 bp sliding windows across each region
bedtools makewindows -b /groups/wyattgrp/projects/chordoma/panel_design/candidates/all_muts_200bp_regions.bed -w 120 -s 1 -i src > all_sliding_windows.bed

# filter for windows that are not edge-clipped (actually 120 bp)
awk '($3-$2) == 120' all_sliding_windows.bed > all_cand_probes.bed

# generate fasta file of from all candidate probes bed file
bedtools getfasta -fi /groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa -name+ -bed all_cand_probes.bed -fo all_cand_seqs.fa 

# run blat with minScore=40 to check for secondary hits with high score in genome
blat /groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa all_cand_seqs.fa -stepSize=5 -repMatch=2253 -minScore=40 -minIdentity=0 blat_all_mins40.psl

# run blat with minScore=30 for detailed filtering of additional criterias
blat /groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa all_cand_seqs.fa -stepSize=5 -repMatch=2253 -minScore=30 -minIdentity=0 blat_all_mins30.psl

# get list of ALL candidate probe_ids
awk 'NR > 5 && ($1 == 120 && $18 == 1)' ./reference/blat_all_mins40.psl | cut -f 10 -d $'\t' > ./reference/all_probe_ids.txt

# get secondary hits that are pslScore > 40
awk '!($1 == 120 && $18 == 1)' ./reference/blat_all_mins40.psl > ./filtering/secondary40_invalid-probes.psl

# get invalid probes that violate: ≤3 contiguous perfect matches (≥20 bp each), and no block >50 bp with >90% identity
python /groups/wyattgrp/users/zshong/codebook//filter_psl.py \
    --path_input ./reference/blat_all_mins30.psl \
    --path_output ./filtering/addbloks_invalid-probes.psl

# get gc content of all probes
seqkit fx2tab ./reference/all_cand_seqs.fa -g -H -n > ./reference/gc_content.tsv

# filter for 30 < gc% < 75 
awk '(30 < $2 && $2 < 75)' ./reference/gc_content.tsv | cut -f 1 -d $'\t' > ./filtering/filtered_gc_ids.txt

# filter for homopolymers >5 bp
python /groups/wyattgrp/users/zshong/codebook/probe_selection/find_homopolymers.py \
    --path_output /groups/wyattgrp/projects/chordoma/panel_design/candidates/filtering/homopolymers_invalid-ids.txt \
    --path_input /groups/wyattgrp/projects/chordoma/panel_design/candidates/reference/all_cand_seqs.fa

# apply no score>40 off-target hits filter
cut -f 10 -d $'\t' ./filtering/secondary40_invalid-probes.psl > ./filtering/off40_invalid_ids.txt
awk 'NR==FNR{a[$1]; next} !($1 in a)' ./filtering/off40_invalid_ids.txt ./reference/all_probe_ids.txt > ./filtering/off40_filtered.txt

# apply additional criteria filter
cut -f 10 -d $'\t' ./filtering/addbloks_invalid-probes.psl > ./filtering/addcritera_invalid_ids.txt
awk 'NR==FNR{a[$1]; next} !($1 in a)' ./filtering/addcritera_invalid_ids.txt ./filtering/off40_filtered.txt > ./filtering/addcritera_off40_filtered.txt

# apply gc content filter
awk 'NR==FNR{a[$1]; next} ($1 in a)' ./filtering/filtered_gc_ids.txt ./filtering/addcritera_off40_filtered.txt > ./qualified/align_crit_passed.txt

# apply homopolymer >5 bp filter
awk 'NR==FNR{a[$1]; next} !($1 in a)' ./filtering/homopolymers_invalid-ids.txt ./qualified/align_crit_passed.txt > ./qualified/final_probes.txt

# for each mutation, get probe in which position of mut is closest to center
python /groups/wyattgrp/users/zshong/codebook/probe_selection/closest_center_probes.py \
    --path_output /groups/wyattgrp/projects/chordoma/panel_design/candidates/qualified \
    --path_input /groups/wyattgrp/projects/chordoma/panel_design/candidates/qualified/final_probes.txt \
    --size 120 # optional, default is 120bp probes

## get probe sequence from probe_id
grep "probe-id" -A 1 /groups/wyattgrp/projects/chordoma/panel_design/candidates/reference/all_cand_seqs.fa
