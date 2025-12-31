#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="/groups/wyattgrp/users/zshong/codebook/probe_selection"
REF_DIR="/groups/wyattgrp/users/jbacon/reference/matti_hg38"
BASE_OUT_DIR="/groups/wyattgrp/projects/chordoma/panel_design"

# Function to process a single TSV
process_tsv() {
    tsv_file="$1"
    base=$(basename "$tsv_file" .tsv)
    OUT_DIR="${BASE_OUT_DIR}/${base}"
    
    echo "Processing $base ..."

    # Create per-file directories
    mkdir -p "${OUT_DIR}/"
    mkdir -p "${OUT_DIR}/filtering"
    mkdir -p "${OUT_DIR}/qualified"

    # Generate 200 bp regions around each mutation
    python "${CODE_DIR}/bed_from_mut_pos.py" \
        --path_output "${OUT_DIR}/" \
        --path_tsv "$tsv_file" \
        --chrom_col "CHROM" \
        --pos_col "POSITION" \
        --id_col "MUTATION ID" \
        --size 100

    # Make 120 bp sliding windows
    bedtools makewindows -b "${OUT_DIR}/all_muts_200bp_regions.bed" -w 120 -s 1 -i src > "${OUT_DIR}/all_sliding_windows.bed"

    # Filter exact 120 bp windows
    awk '($3-$2) == 120' "${OUT_DIR}/all_sliding_windows.bed" > "${OUT_DIR}/all_cand_probes.bed"

    # Extract sequences
    bedtools getfasta -fi "${REF_DIR}/hg38_masked.fa" -name+ -bed "${OUT_DIR}/all_cand_probes.bed" -fo "${OUT_DIR}/all_cand_seqs.fa"

    # Run BLAT
    blat "${REF_DIR}/hg38_masked.fa" "${OUT_DIR}/all_cand_seqs.fa" -stepSize=5 -repMatch=2253 -minScore=40 -minIdentity=0 "${OUT_DIR}/blat_all_mins40.psl"

    # Extract valid probe IDs
    awk 'NR>5 && ($1==120 && $18==1)' "${OUT_DIR}/blat_all_mins40.psl" | cut -f10 > "${OUT_DIR}/all_probe_ids.txt"

    # Secondary hits
    awk '!($1==120 && $18==1)' "${OUT_DIR}/blat_all_mins40.psl" > "${OUT_DIR}/filtering/secondary40_invalid-probes.psl"

    # Additional criteria filter
    python "${CODE_DIR}/filter_psl.py" \
        --path_input "${OUT_DIR}/blat_all_mins40.psl" \
        --path_output "${OUT_DIR}/filtering/addbloks_invalid-probes.psl"

    # GC content filter
    seqkit fx2tab "${OUT_DIR}/all_cand_seqs.fa" -g -H -n > "${OUT_DIR}/gc_content.tsv"
    awk '(30 < $2 && $2 < 75)' "${OUT_DIR}/gc_content.tsv" | cut -f1 > "${OUT_DIR}/filtering/filtered_gc_ids.txt"

    # Homopolymer filter
    python "${CODE_DIR}/find_homopolymers.py" \
        --path_output "${OUT_DIR}/filtering/homopolymers_invalid-ids.txt" \
        --path_input "${OUT_DIR}/all_cand_seqs.fa"

    # Apply off-target filter
    cut -f10 "${OUT_DIR}/filtering/secondary40_invalid-probes.psl" > "${OUT_DIR}/filtering/off40_invalid_ids.txt"
    awk 'NR==FNR{a[$1]; next} !($1 in a)' "${OUT_DIR}/filtering/off40_invalid_ids.txt" "${OUT_DIR}/all_probe_ids.txt" > "${OUT_DIR}/filtering/off40_filtered.txt"

    # Apply additional criteria
    cut -f10 "${OUT_DIR}/filtering/addbloks_invalid-probes.psl" > "${OUT_DIR}/filtering/addcritera_invalid_ids.txt"
    awk 'NR==FNR{a[$1]; next} !($1 in a)' "${OUT_DIR}/filtering/addcritera_invalid_ids.txt" "${OUT_DIR}/filtering/off40_filtered.txt" > "${OUT_DIR}/filtering/addcritera_off40_filtered.txt"

    # Apply GC filter
    awk 'NR==FNR{a[$1]; next} ($1 in a)' "${OUT_DIR}/filtering/filtered_gc_ids.txt" "${OUT_DIR}/filtering/addcritera_off40_filtered.txt" > "${OUT_DIR}/qualified/align_crit_passed.txt"

    # Apply homopolymer filter
    awk 'NR==FNR{a[$1]; next} !($1 in a)' "${OUT_DIR}/filtering/homopolymers_invalid-ids.txt" "${OUT_DIR}/qualified/align_crit_passed.txt" > "${OUT_DIR}/qualified/final_probes.txt"

    # Closest center probes
    python "${CODE_DIR}/closest_center_probes.py" \
        --path_output "${OUT_DIR}/qualified" \
        --path_input "${OUT_DIR}/qualified/final_probes.txt" \
        --size 120

    # Merge probe_id back into original mutation TSV
    python - <<END
import pandas as pd
muts = pd.read_csv("$tsv_file", sep="\t")
probes = pd.read_csv("${OUT_DIR}/qualified/probe_loc_stats.tsv", sep="\t")
muts = muts.merge(probes[['Sample','Chrom','Position','probe_id']], on=['Sample','Chrom','Position'], how='left')
muts.to_csv("${OUT_DIR}/qualified/${base}_with_probes.tsv", sep="\t", index=False)
END

    echo "Finished processing $base"
}

export -f process_tsv
export CODE_DIR REF_DIR BASE_OUT_DIR

# Run on all TSV files in parallel
process_tsv /groups/wyattgrp/projects/chordoma/mutations/TWIST_std_muts.tsv