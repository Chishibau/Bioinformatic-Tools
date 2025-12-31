#!/usr/bin/env nextflow

/*
 *  PROBE SELECTION PIPELINE — NEXTFLOW VERSION
 *  Includes:
 *    - Generate 200 bp regions around mutations
 *    - Sliding windows (120 bp)
 *    - BLAT filtering (minScore=40)
 *    - Off-target filtering
 *    - GC-content filtering
 *    - Homopolymer filtering
 *    - Closest-center probe selection
 *    - Merge probe_id into original mutations file
 */

params.mutations   = "/groups/wyattgrp/projects/brachy_tracks/mutations/mutato_analyze/muts_somatic_mutato_analyze.tsv"
params.ref_fasta   = "/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
params.output      = "./results"
params.size        = 100
params.sliding     = 120

process MAKE_REGIONS {
    publishDir "${params.output}/regions", mode: 'copy'

    input:
    path muts_file = params.mutations

    output:
    path "all_muts_200bp_regions.bed"

    """
    python /groups/wyattgrp/users/zshong/codebook/probe_selection/bed_from_mut_pos.py \
        --path_output . \
        --path_tsv ${muts_file} \
        --chrom_col "Chrom" \
        --pos_col "Position" \
        --id_col "Sample" \
        --size ${params.size}
    """
}

process SLIDING_WINDOWS {
    publishDir "${params.output}/windows", mode: 'copy'

    input:
    path bed_regions

    output:
    path "all_cand_probes.bed"

    """
    bedtools makewindows \
        -b ${bed_regions} \
        -w ${params.sliding} \
        -s 1 \
        -i src > all_sliding_windows.bed

    awk '(\$3-\$2)==${params.sliding}' all_sliding_windows.bed > all_cand_probes.bed
    """
}

process GET_FASTA {
    publishDir "${params.output}/fasta", mode: 'copy'

    input:
    path bed = "all_cand_probes.bed"

    output:
    path "all_cand_seqs.fa"

    """
    bedtools getfasta \
        -fi ${params.ref_fasta} \
        -name+ \
        -bed ${bed} \
        -fo all_cand_seqs.fa
    """
}

process RUN_BLAT {
    publishDir "${params.output}/blat", mode: 'copy'

    input:
    path fasta = "all_cand_seqs.fa"

    output:
    path "blat_all_mins40.psl"

    """
    blat ${params.ref_fasta} ${fasta} \
        -stepSize=5 -repMatch=2253 -minScore=40 -minIdentity=0 \
        blat_all_mins40.psl
    """
}

process FILTER_BLAT {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    path psl = "blat_all_mins40.psl"

    output:
    path "all_probe_ids.txt"
    path "secondary40_invalid-probes.psl"

    """
    # valid probes
    awk 'NR > 5 && (\$1 == 120 && \$18 == 1)' ${psl} | cut -f10 > all_probe_ids.txt

    # invalid: secondary hits score>40
    awk '! (\$1 == 120 && \$18 == 1)' ${psl} > secondary40_invalid-probes.psl
    """
}

process ADDITIONAL_FILTERS {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    path psl = "blat_all_mins40.psl"

    output:
    path "addcritera_invalid_ids.txt"

    """
    python /groups/wyattgrp/users/zshong/codebook/probe_selection/filter_psl.py \
        --path_input ${psl} \
        --path_output addbloks_invalid-probes.psl

    cut -f 10 addbloks_invalid-probes.psl > addcritera_invalid_ids.txt
    """
}

process GC_CONTENT {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    path fasta = "all_cand_seqs.fa"

    output:
    path "filtered_gc_ids.txt"

    """
    seqkit fx2tab ${fasta} -g -H -n > gc_content.tsv
    awk '(30 < \$2 && \$2 < 75) {print \$1}' gc_content.tsv \
        > filtered_gc_ids.txt
    """
}

process HOMOPOLYMERS {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    path fasta = "all_cand_seqs.fa"

    output:
    path "homopolymers_invalid-ids.txt"

    """
    python /groups/wyattgrp/users/zshong/codebook/probe_selection/find_homopolymers.py \
        --path_output homopolymers_invalid-ids.txt \
        --path_input ${fasta}
    """
}

process APPLY_FILTERS {
    publishDir "${params.output}/qualified", mode: 'copy'

    input:
    path all_ids
    path invalid40
    path invalid_crit
    path valid_gc
    path invalid_hpoly

    output:
    path "final_probes.txt"

    """
    # remove >40 score off-targets
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' ${invalid40} ${all_ids} > off40_filtered.txt

    # remove add. criteria-invalid
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' ${invalid_crit} off40_filtered.txt > addcritera_off40_filtered.txt

    # keep GC good
    awk 'NR==FNR{good[\$1];next} (\$1 in good)' ${valid_gc} addcritera_off40_filtered.txt > gc_passed.txt

    # remove homopolymers invalid
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' ${invalid_hpoly} gc_passed.txt > final_probes.txt
    """
}

process CLOSEST_CENTER {
    publishDir "${params.output}/qualified", mode: 'copy'

    input:
    path final_ids = "final_probes.txt"

    output:
    path "probe_loc_stats.tsv"

    """
    python /groups/wyattgrp/users/zshong/codebook/probe_selection/closest_center_probes.py \
        --path_output . \
        --path_input ${final_ids} \
        --size ${params.sliding}
    """
}

process MERGE_PROBES {
    publishDir "${params.output}/final", mode: 'copy'

    input:
    path locstats = "probe_loc_stats.tsv"
    path muts = params.mutations

    output:
    path "mutations_with_probe_id.tsv"

    """
    python - << 'EOF'
import pandas as pd

muts = pd.read_csv("${muts}", sep="\\t")
probes = pd.read_csv("${locstats}", sep="\\t")

merged = muts.merge(
    probes[['Sample','Chrom','Position','probe_id']],
    on=['Sample','Chrom','Position'],
    how='left'
)

merged.to_csv("mutations_with_probe_id.tsv", sep="\\t", index=False)
EOF
    """
}

workflow {
    Channel
        .fromPath(params.mutations)
        | MAKE_REGIONS
        | SLIDING_WINDOWS
        | GET_FASTA
        | RUN_BLAT

    RUN_BLAT.out | FILTER_BLAT
    RUN_BLAT.out | ADDITIONAL_FILTERS
    GET_FASTA.out | GC_CONTENT
    GET_FASTA.out | HOMOPOLYMERS

    APPLY_FILTERS(
        FILTER_BLAT.out[0],
        FILTER_BLAT.out[1],
        ADDITIONAL_FILTERS.out,
        GC_CONTENT.out,
        HOMOPOLYMERS.out
    ) | CLOSEST_CENTER | MERGE_PROBES
}
