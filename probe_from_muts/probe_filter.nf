#!/usr/bin/env nextflow

/*
 *  PROBE SELECTION PIPELINE — NEXTFLOW VERSION
 *  Includes:
 *    - Generate defined bp regions around mutations
 *    - Sliding windows (specified probe length in bp)
 *    - BLAT filtering (minScore for secondary alignment provided)
 *    - Off-target filtering
 *    - GC-content filtering
 *    - Homopolymer filtering
 *    - Closest-center probe selection
 *    - Merge probe_id into original mutations file
 *    - Outputs mutations that are closely positioned (within half of probe-length)
 */

params.mutations   = null
params.ref_fasta   = "/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
params.output      = "./probe_filtering"
params.size        = 100
params.probe_size  = 120
params.blat_threshold = 40

process MAKE_REGIONS {
    publishDir "${params.output}/inter_files", mode: 'copy'

    input:
    path muts_file = params.mutations

    output:
    path "all_muts_200bp_regions.bed"

    """
    bed_from_mut_pos.py \
        --path_output . \
        --path_tsv ${muts_file} \
        --chrom_col "Chrom" \
        --pos_col "Position" \
        --id_col "Sample" \
        --size ${params.size}
    """
}

process SLIDING_WINDOWS {
    publishDir "${params.output}/inter_files", mode: 'copy'

    input:
    path bed_regions

    output:
    path "all_cand_probes.bed"

    """
    bedtools makewindows \
        -b ${bed_regions} \
        -w ${params.probe_size} \
        -s 1 \
        -i src > all_sliding_windows.bed

    awk '(\$3-\$2)==${params.probe_size}' all_sliding_windows.bed > all_cand_probes.bed
    """
}

process GET_FASTA {
    publishDir "${params.output}/inter_files", mode: 'copy'

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
    publishDir "${params.output}/inter_files", mode: 'copy'

    input:
    path fasta = "all_cand_seqs.fa"

    output:
    path "blat_all_mins${blat_threshold}.psl"

    """
    blat ${params.ref_fasta} ${fasta} \
        -stepSize=5 -repMatch=2253 -minScore=${blat_threshold} -minIdentity=0 \
        blat_all_mins${blat_threshold}.psl
    """
}

process FILTER_BLAT {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    path psl = "blat_all_mins${blat_threshold}.psl"

    output:
    path "all_probe_ids.txt"
    path "secondary${blat_threshold}_invalid-probes.psl"

    """
    # valid probes
    awk 'NR > 5 && (\$1 == 120 && \$18 == 1)' ${psl} | cut -f10 > all_probe_ids.txt

    # invalid: secondary hits score>40
    awk '! (\$1 == 120 && \$18 == 1)' ${psl} > secondary${blat_threshold}_invalid-probes.psl
    """
}

process ADDITIONAL_FILTERS {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    path psl = "blat_all_mins${blat_threshold}.psl"

    output:
    path "addcritera_invalid_ids.txt"

    """
    filter_psl.py \
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
    find_homopolymers.py \
        --path_output homopolymers_invalid-ids.txt \
        --path_input ${fasta}
    """
}

process APPLY_FILTERS {
    publishDir "${params.output}/qualified", mode: 'copy'

    input:
    path all_ids
    path invalid_secondary
    path invalid_crit
    path valid_gc
    path invalid_hpoly

    output:
    path "final_probes.txt"

    """
    # remove >40 score off-targets
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' ${invalid_secondary} ${all_ids} > no_off${blat_threshold}.txt

    # remove muts with perfect match secondary alignments
    sort final_probes.txt | uniq -u > no_secondary.txt

    # remove add. criteria-invalid
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' ${invalid_crit} no_secondary.txt > addcritera_secondary_filtered.txt

    # keep GC good
    awk 'NR==FNR{good[\$1];next} (\$1 in good)' ${valid_gc} addcritera_secondary_filtered.txt > gc_passed.txt

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
    closest_center_probes.py \
        --path_output . \
        --path_input ${final_ids} \
        --size ${params.probe_size}
    """
}

process MERGE_PROBES {
    publishDir "${params.output}/results", mode: 'copy'

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

        merged.to_csv("muts_with_probe_label.tsv", sep="\\t", index=False)
        EOF
    """
}

process CHECK_CLOSE_MUTS {
    publishDir "${params.output}/results", mode: 'copy'

    input:
    path locstats = "probe_loc_stats.tsv"

    output:
    path "close_loc_muts.tsv"
    path "no_close_muts.tsv"

    """
    check_close_muts.py \
        --path_output . \
        --path_input ${locstats} \
        --size ${params.probe_size}
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
    ) | CLOSEST_CENTER | MERGE_PROBES | CHECK_CLOSE_MUTS
}
