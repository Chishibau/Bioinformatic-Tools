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

params.regions   = null
params.ref_fasta   = "/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
params.output      = "./region_tiling"
params.probe_size  = 120
params.blat_threshold = 40
params.bp_homopoly = 5
params.gc_upper = 75
params.gc_lower = 30

process SLIDING_WINDOWS {
    publishDir "${params.output}/inter_files", mode: 'copy'

    input:
    tuple val(chr), val(start), val(end), val(id), val(sliding), val(spacing)

    output:
    tuple val(id), val(spacing), path("${id}_cand_probes.bed")

    script:
    """
    cat << EOF > region.bed
    ${chr}\t${start}\t${end}\t${id}
    EOF

    bedtools makewindows \
        -b region.bed  \
        -w ${params.probe_size} \
        -s ${sliding} \
        -i src > ${id}_windows.bed

    awk '(\$3-\$2)==${params.probe_size}' ${id}_windows.bed > ${id}_cand_probes.bed
    """
}

process GET_FASTA {
    publishDir "${params.output}/inter_files", mode: 'copy'

    input:
    tuple val(id), val(spacing), path(bed)

    output:
    tuple val(id), val(spacing), path("${id}_cand_seqs.fa")

    script:
    """
    bedtools getfasta \
        -fi ${params.ref_fasta} \
        -name+ \
        -bed ${bed} \
        -fo ${id}_cand_seqs.fa
    """
}

process RUN_BLAT {
    publishDir "${params.output}/inter_files", mode: 'copy'

    input:
    tuple val(id), val(spacing), path(fasta)

    output:
    tuple val(id), val(spacing), path("${id}_blat_mins${params.blat_threshold}.psl")

    script:
    """
    blat ${params.ref_fasta} ${fasta} \
        -stepSize=5 -repMatch=2253 -minScore=${params.blat_threshold} -minIdentity=0 \
        ${id}_blat_mins${params.blat_threshold}.psl
    """
}

process FILTER_BLAT {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    tuple val(id), val(spacing), path(psl)

    output:
    tuple val(id), val(spacing), path("${id}_all_ids.txt"), path("${id}_blat${params.blat_threshold}_invalids.txt")

    script:
    """
    # Valid: perfect score (==probe_size) with exactly 1 alignment block
    awk 'NR > 5 && (\$1 == ${params.probe_size} && \$18 == 1)' ${psl} \
        | cut -f 10 \
        | sort -u \
        > ${id}_all_ids.txt

    # Invalid: any alignment that is NOT a perfect single-block hit
    awk 'NR > 5 && !(\$1 == ${params.probe_size} && \$18 == 1)' ${psl} \
        | cut -f 10 \
        | sort -u \
        > ${id}_blat${params.blat_threshold}_invalids.txt
    """
}

process ADDITIONAL_FILTERS {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    tuple val(id), val(spacing), path(psl)

    output:
    tuple val(id), val(spacing), path("${id}_addcrit_invalids.txt")

    script:
    """
    filter_psl.py \
        --path_input ${psl} \
        --path_output ${id}_addbloks_invalid_probes.psl

    cut -f 10 ${id}_addbloks_invalid_probes.psl > ${id}_addcrit_invalids.txt
    """
}

process GC_CONTENT {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    tuple val(id), val(spacing), path(fasta)

    output:
    tuple val(id), val(spacing), path("${id}_qual_gc_ids.txt")

    script:
    """
    seqkit fx2tab ${fasta} -g -H -n > ${id}_gc_content.tsv
    awk '(${params.gc_lower} < \$2 && \$2 < ${params.gc_upper}) {print \$1}' ${id}_gc_content.tsv \
        > ${id}_qual_gc_ids.txt
    """
}

process HOMOPOLYMERS {
    publishDir "${params.output}/filtering", mode: 'copy'

    input:
    tuple val(id), val(spacing), path(fasta)

    output:
    tuple val(id), val(spacing), path("${id}_homopoly_invalids.txt")

    script:
    """
    find_homopolymers.py \
        --path_output ${id}_homopoly_invalids.txt \
        --path_input ${fasta} \
        --bp ${params.bp_homopoly}
    """
}

process APPLY_FILTERS {
    publishDir "${params.output}/qualified", mode: 'copy'

    input:
    tuple val(id), val(spacing), path(all_ids), path(invalid_secondary_ids), path(invalid_crit), path(valid_gc), path(invalid_hpoly)

    output:
    tuple val(id), val(spacing), path("${id}_qual_probes.txt")

    script:
    """
    # Remove probes with off-target alignments above score threshold
    if [ -s ${invalid_secondary_ids} ]; then
    awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' \
        ${invalid_secondary_ids} ${all_ids} \
        > ${id}_no_off${params.blat_threshold}.txt
    else 
        cp ${all_ids} ${id}_no_off${params.blat_threshold}.txt
    fi

    # Remove probes failing additional criteria
    if [ -s ${invalid_crit} ]; then
        awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' \
            ${invalid_crit} ${id}_no_off${params.blat_threshold}.txt \
            > ${id}_addcriteria_filtered.txt
    else
        cp ${id}_no_off${params.blat_threshold}.txt ${id}_addcriteria_filtered.txt
    fi

    # Keep only probes passing GC filter
    if [ -s ${valid_gc} ]; then
        awk 'NR==FNR{good[\$1];next} (\$1 in good)' \
            ${valid_gc} ${id}_addcriteria_filtered.txt \
            > ${id}_gc_passed.txt
    else
        cp ${id}_addcriteria_filtered.txt ${id}_gc_passed.txt
    fi

    # Remove probes with homopolymer runs
    if [ -s ${invalid_hpoly} ]; then
        awk 'NR==FNR{bad[\$1];next} !(\$1 in bad)' \
            ${invalid_hpoly} ${id}_gc_passed.txt \
            > ${id}_qual_probes.txt
    else
        cp ${id}_gc_passed.txt ${id}_qual_probes.txt
    fi
    """
}

process TILE_REGION {
    publishDir "${params.output}/results", mode: 'copy'

    input:
    tuple val(id), val(spacing), path(final_ids)

    output:
    path("${id}_probes.bed")

    script:
    """
    tile_region.py \
        --output ${id}_probes.bed \
        --input ${final_ids} \
        --bp ${spacing}
    """
}

workflow {
    // create input channels for each row/region in tsv file
    if (params.regions) {
        regions_ch = Channel.fromPath(params.regions)
            .splitCsv(header: true, sep: '\t')
            .map { row -> 
                tuple(
                    row.chr, 
                    row.start, 
                    row.end, 
                    row.id, 
                    row.sliding, 
                    row.spacing
                ) 
            }
    } else {
        error "Please provide --regions (single file)"
    }

    // generate candidate probes and blat sequences
    fasta_ch = regions_ch | SLIDING_WINDOWS | GET_FASTA
    blat_ch = RUN_BLAT(fasta_ch)
    
    // perform filtering
    filter_blat_ch = FILTER_BLAT(blat_ch)
    add_filter_ch = ADDITIONAL_FILTERS(blat_ch)
    gc_ch = GC_CONTENT(fasta_ch)
    hpoly_ch = HOMOPOLYMERS(fasta_ch)
    
    // combine filtered ids into 1 channel
    combined_ch = filter_blat_ch       // (sample_id, all_ids, invalid_secondary_ids)
    .join(add_filter_ch, by: [0, 1])               // + invalid_crit
    .join(gc_ch, by: [0, 1])                       // + valid_gc
    .join(hpoly_ch, by: [0, 1])                    // + invalid_hpoly
    
    stats_ch = combined_ch | APPLY_FILTERS | TILE_REGION
}