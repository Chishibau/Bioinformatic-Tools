#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.mutations_tsv = null
params.bam_dir       = null
params.reference     = "/groups/wyattgrp/reference/hg38/hg38.fa"
params.outdir        = "results"

if (!params.mutations_tsv) error "Provide --mutations_tsv"
if (!params.bam_dir) error "Provide --bam_dir"

process PARSE_MUTATIONS {
    publishDir "${params.outdir}/split_mutations", mode: 'copy'

    input:
    path mutations_tsv

    output:
    path "patient_*.tsv"

    script:
    """
    python3 - << 'EOF'
    import pandas as pd

    df = pd.read_csv("${mutations_tsv}", sep='\\t')
    df = df.rename(columns={'ID':'patient_id','ref':'ref_allele','alt':'alt_allele'})

    if 'alt_allele' not in df.columns:
        df['alt_allele'] = ''
    if 'origin' not in df.columns:
        df['origin'] = 'unknown'

    patient_df = df[df['origin']=='patient']
    shared_df  = df[df['origin']!='patient']

    for pid in patient_df['patient_id'].dropna().unique():
        safe = str(pid).replace('/', '_')
        subset = patient_df[patient_df['patient_id']==pid]
        combined = pd.concat([subset, shared_df], ignore_index=True)
        combined = combined.drop_duplicates(subset=['chrom','position','ref_allele','alt_allele'])
        combined.to_csv(f"patient_{safe}.tsv", sep='\\t', index=False)
    EOF
    """
}

process TSV_TO_BED {
    tag "${patient_id}"

    input:
    tuple val(patient_id), path(tsv)

    output:
    tuple val(patient_id), path(tsv), path("regions.bed")

    script:
    """
    python3 - << 'EOF'
    import pandas as pd

    df = pd.read_csv("${tsv}", sep='\\t')

    if df.empty:
        open("regions.bed","w").close()
    else:
        bed = pd.DataFrame({
            'chrom': df['chrom'],
            'start': df['position'] - 1,
            'end': df['position']
        }).drop_duplicates()
        
        bed.to_csv("regions.bed", sep='\\t', header=False, index=False)
    EOF
    """
}

process FIND_BAM {
    tag "${patient_id}"
    
    input:
    tuple val(patient_id), path(tsv), path(bed)
    val bam_dir
    
    output:
    tuple val(patient_id), path(tsv), path(bed), env(BAM_PATH), env(BAM_INDEX_PATH)
    
    script:
    """
    #!/bin/bash
    
    # Find BAM file matching patient ID
    BAM_PATH=\$(find ${bam_dir} -name "*${patient_id}*FFPE*.bam" | head -1)
    
    if [ -z "\$BAM_PATH" ]; then
        echo "ERROR: No BAM file found for patient ${patient_id}" >&2
        exit 1
    fi
    
    if [ ! -f "\$BAM_PATH" ]; then
        echo "ERROR: BAM file does not exist: \$BAM_PATH" >&2
        exit 1
    fi
    
    echo "Found BAM for ${patient_id}: \$BAM_PATH"
    
    # Find BAM index
    if [ -f "\${BAM_PATH}.bai" ]; then
        BAM_INDEX_PATH="\${BAM_PATH}.bai"
    elif [ -f "\${BAM_PATH%.*}.bai" ]; then
        BAM_INDEX_PATH="\${BAM_PATH%.*}.bai"
    else
        echo "ERROR: No BAM index found for \$BAM_PATH" >&2
        echo "Please index the BAM file with: samtools index \$BAM_PATH" >&2
        exit 1
    fi
    
    echo "Found BAM index: \$BAM_INDEX_PATH"
    
    # Export for Nextflow to capture
    echo \$BAM_PATH
    echo \$BAM_INDEX_PATH
    """
}

process BCFTOOLS_VAF {
    tag "${patient_id}"
    publishDir "${params.outdir}/vaf", mode: 'copy'

    input:
    tuple val(patient_id), path(tsv), path(bed), val(bam_path), val(bam_index), path(reference)

    output:
    path "${patient_id}_vaf.tsv"

    script:
    """
    bcftools mpileup \
        -f ${reference} \
        -R ${bed} \
        -B \
        -d 1600 \
        -a FORMAT/AD,FORMAT/DP \
        -Ou ${bam_path} | \
    bcftools query \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP]\\t[%AD]\\n' \
        > raw_counts.tsv

    python3 - << 'EOF'
    import pandas as pd

    mut_df = pd.read_csv("${tsv}", sep='\\t')
    bam_path = "${bam_path}"

    cols = ['chrom','position','ref','alt_list','depth','ad']
    try:
        df = pd.read_csv("raw_counts.tsv", sep='\\t', names=cols)
    except:
        df = pd.DataFrame(columns=cols)

    rows = []
    for _, r in df.iterrows():
        alts = str(r.alt_list).split(',') if pd.notna(r.alt_list) else []
        counts = [int(x) if x.isdigit() else 0 for x in str(r.ad).split(',')] if pd.notna(r.ad) else []
        
        rows.append({
            'chrom': r.chrom,
            'position': int(r.position),
            'ref': str(r.ref).upper(),
            'alts': alts,
            'ref_count': counts[0] if counts else 0,
            'alt_counts': counts[1:] if len(counts)>1 else [],
            'depth': sum(counts)
        })

    bcf_df = pd.DataFrame(rows)

    final = []
    for _, mut in mut_df.iterrows():
        chrom, pos, origin = str(mut.chrom), int(mut.position), str(mut.origin)
        alt = str(mut.alt_allele).upper() if pd.notna(mut.alt_allele) else ""
        
        match = bcf_df[(bcf_df.chrom==chrom)&(bcf_df.position==pos)]
        
        if match.empty:
            final.append(["${patient_id}", chrom, pos, 0, 0, 0, alt, 0.0, origin, bam_path])
            continue
        
        row = match.iloc[0]
        alt_count = row.alt_counts[row.alts.index(alt)] if alt in row.alts else 0
        depth = row.depth
        
        final.append(["${patient_id}", chrom, pos, depth, row.ref_count, alt_count, alt, round((alt_count/depth*100) if depth>0 else 0, 2), origin, bam_path])

    out = pd.DataFrame(final, columns=['patient_id', 'chrom','position','depth','ref_count','alt_count','alt','alt_VAF', 'origin', 'bam_file'])
    out.to_csv("${patient_id}_vaf.tsv", sep='\\t', index=False)
    EOF
    """
}

process COMBINE_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path vafs

    output:
    path "final_vaf_results.tsv"

    script:
    """
    python3 - << 'EOF'
    import pandas as pd
    import glob

    files = glob.glob("*_vaf.tsv")
    dfs = [pd.read_csv(f, sep='\\t') for f in files]

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        combined = combined.sort_values(['patient_id', 'chrom', 'position'])
        combined.to_csv("final_vaf_results.tsv", sep='\\t', index=False)
    else:
        # Create empty file with header
        pd.DataFrame(columns=['patient_id','bam_file','chrom','position','depth','ref_count','alt_count','alt','alt_VAF','origin']).to_csv("final_vaf_results.tsv", sep='\\t', index=False)
    EOF
    """
}

workflow {
    // Input files
    reference = file(params.reference)
    bam_dir = params.bam_dir
    
    // Step 1: Parse mutations and split by patient
    patient_tsvs = PARSE_MUTATIONS(params.mutations_tsv)
        .flatten()
        .map { tsv -> 
            def patient_id = tsv.baseName.replace("patient_", "")
            tuple(patient_id, tsv)
        }
    
    // Step 2: Convert TSV to BED
    patient_beds = TSV_TO_BED(patient_tsvs)
    
    // Step 3: Find matching BAM file for each patient
    patient_bams = FIND_BAM(patient_beds, bam_dir)
    
    // Step 4: Calculate VAF
    vaf_input = patient_bams.map { patient_id, tsv, bed, bam_path, bam_index ->
        tuple(patient_id, tsv, bed, bam_path, bam_index, reference)
    }
    
    vaf_results = BCFTOOLS_VAF(vaf_input)
    
    // Step 5: Combine results
    COMBINE_RESULTS(vaf_results.collect())
}

workflow.onComplete {
    log.info """
        Pipeline completed!
        Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Results: ${params.outdir}/final_vaf_results.tsv
        """
        .stripIndent()
}