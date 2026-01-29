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
params.mutations_dir = null
params.ref_fasta   = "/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
params.output      = "./probe_filtering"
params.size        = 100
params.probe_size  = 120
params.blat_threshold = 40

process MAKE_REGIONS {
    tag "${muts_file.simpleName}"
    publishDir "${params.output}/${muts_file.simpleName}/inter_files", mode: 'copy'

    input:
    path muts_file

    output:
    tuple val("${muts_file.simpleName}"), path("all_muts_200bp_regions.bed")

    script:
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
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/inter_files", mode: 'copy'

    input:
    tuple val(sample_id), path(bed_regions)

    output:
    tuple val(sample_id), path("all_cand_probes.bed")

    script:
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
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/inter_files", mode: 'copy'

    input:
    tuple val(sample_id), path(bed)

    output:
    tuple val(sample_id), path("all_cand_seqs.fa")

    script:
    """
    bedtools getfasta \
        -fi ${params.ref_fasta} \
        -name+ \
        -bed ${bed} \
        -fo all_cand_seqs.fa
    """
}

process RUN_BLAT {
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/inter_files", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("blat_all_mins${params.blat_threshold}.psl")

    script:
    """
    blat ${params.ref_fasta} ${fasta} \
        -stepSize=5 -repMatch=2253 -minScore=${params.blat_threshold} -minIdentity=0 \
        blat_all_mins${params.blat_threshold}.psl
    """
}

process FILTER_BLAT {
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/filtering", mode: 'copy'

    input:
    tuple val(sample_id), path(psl)

    output:
    tuple val(sample_id), path("all_probe_ids.txt"), path("secondary${params.blat_threshold}_invalid-probes.psl")

    script:
    """
    # valid probes
    awk 'NR > 5 && (\$1 == 120 && \$18 == 1)' ${psl} | cut -f 10 > all_probe_ids.txt

    # invalid: secondary hits score>40
    awk '! (\$1 == 120 && \$18 == 1)' ${psl} > secondary${params.blat_threshold}_invalid-probes.psl
    """
}

process ADDITIONAL_FILTERS {
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/filtering", mode: 'copy'

    input:
    tuple val(sample_id), path(psl)

    output:
    tuple val(sample_id), path("addcritera_invalid_ids.txt")

    script:
    """
    filter_psl.py \
        --path_input ${psl} \
        --path_output addbloks_invalid-probes.psl

    cut -f 10 addbloks_invalid-probes.psl > addcritera_invalid_ids.txt
    """
}

process GC_CONTENT {
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/filtering", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("filtered_gc_ids.txt")

    script:
    """
    seqkit fx2tab ${fasta} -g -H -n > gc_content.tsv
    awk '(30 < \$2 && \$2 < 75) {print \$1}' gc_content.tsv \
        > filtered_gc_ids.txt
    """
}

process HOMOPOLYMERS {
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/filtering", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("homopolymers_invalid-ids.txt")

    script:
    """
    find_homopolymers.py \
        --path_output homopolymers_invalid-ids.txt \
        --path_input ${fasta}
    """
}

process APPLY_FILTERS {
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/qualified", mode: 'copy'

    input:
    tuple val(sample_id), path(all_ids), path(invalid_secondary), path(invalid_crit), path(valid_gc), path(invalid_hpoly)

    output:
    tuple val(sample_id), path("final_probes.txt")

    script:
    """
    # Debug: check input file sizes
    echo "Input file sizes:"
    wc -l ${all_ids} ${invalid_secondary} ${invalid_crit} ${valid_gc} ${invalid_hpoly}
    
    # remove >40 score off-targets
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' ${invalid_secondary} ${all_ids} > no_off${params.blat_threshold}.txt

    # remove muts with perfect match secondary alignments
    sort no_off${params.blat_threshold}.txt | uniq -u > no_secondary.txt

    # remove add. criteria-invalid
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' ${invalid_crit} no_secondary.txt > addcritera_secondary_filtered.txt

    # keep GC good
    awk 'NR==FNR{good[\$1];next} (\$1 in good)' ${valid_gc} addcritera_secondary_filtered.txt > gc_passed.txt

    # remove homopolymers invalid
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' ${invalid_hpoly} gc_passed.txt > final_probes.txt
    
    # If empty, create with at least a comment for debugging
    if [ ! -s final_probes.txt ]; then
        echo "# No probes passed all filters" > final_probes.txt
    fi
    """
}

process CLOSEST_CENTER {
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/qualified", mode: 'copy'

    input:
    tuple val(sample_id), path(final_ids)

    output:
    tuple val(sample_id), path("probe_loc_stats.tsv")

    script:
    """
    closest_center_probes.py \
        --path_output . \
        --path_input ${final_ids} \
        --size ${params.probe_size}
    """
}

process MERGE_PROBES {
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/results", mode: 'copy'

    input:
    tuple val(sample_id), path(locstats), path(muts)

    output:
    tuple val(sample_id), path("muts_with_probe_label.tsv")

    script:
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
    tag "${sample_id}"
    publishDir "${params.output}/${sample_id}/results", mode: 'copy'

    input:
    tuple val(sample_id), path(locstats)

    output:
    tuple val(sample_id), path("close_loc_muts.tsv"), path("no_close_muts.tsv")

    script:
    """
    check_close_muts.py \
        --path_output . \
        --path_input ${locstats} \
        --size ${params.probe_size}
    """
}

workflow {
    // Create input channel from either single file or directory
    if (params.mutations_dir) {
        muts_ch = Channel.fromPath("${params.mutations_dir}/*.{tsv,txt,csv}")
    } else if (params.mutations) {
        muts_ch = Channel.fromPath(params.mutations)
    } else {
        error "Please provide either --mutations (single file) or --mutations_dir (directory)"
    }

    // Main pipeline
    fasta_ch = muts_ch | MAKE_REGIONS | SLIDING_WINDOWS | GET_FASTA
    blat_ch = RUN_BLAT(fasta_ch)
    
    // Get all filter outputs
    filter_blat_ch = FILTER_BLAT(blat_ch)
    add_filter_ch = ADDITIONAL_FILTERS(blat_ch)
    gc_ch = GC_CONTENT(fasta_ch)
    hpoly_ch = HOMOPOLYMERS(fasta_ch)
    
    // Combine all filters by sample_id - keep all 3 elements from filter_blat_ch
    combined_ch = filter_blat_ch
        .join(add_filter_ch)
        .join(gc_ch)
        .join(hpoly_ch)
        .map { sample_id, all_ids, invalid_sec, invalid_crit, valid_gc, invalid_hpoly ->
            tuple(sample_id, all_ids, invalid_sec, invalid_crit, valid_gc, invalid_hpoly)
        }
    
    // Final steps
    stats_ch = combined_ch | APPLY_FILTERS | CLOSEST_CENTER
    
    MERGE_PROBES(stats_ch.join(muts_ch.map { file -> tuple(file.simpleName, file) }))
    CHECK_CLOSE_MUTS(stats_ch)
}