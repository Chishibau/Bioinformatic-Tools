#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.mutations_tsv = null
params.bam_dir       = null
params.reference     = "/groups/wyattgrp/reference/hg38/hg38.fa"
params.outdir        = "results"
params.controls      = false

if( !params.mutations_tsv )
    error "Provide --mutations_tsv"

if( !params.bam_dir )
    error "Provide --bam_dir"


// FILTER MUTATIONS  (strip indels: any row where len(ref) != 1 or len(alt) != 1)
process FILTER_MUTATIONS {

    publishDir "${params.outdir}/mutations", mode: 'copy'

    input:
    path mutations_tsv

    output:
    path "mutations_snv.tsv"

    script:
    """
    python3 - << 'EOF'
    import pandas as pd

    df = pd.read_csv("${mutations_tsv}", sep='\\t')

    before = len(df)

    df = df[
        (df["ref"].astype(str).str.len() == 1) &
        (df["alt"].astype(str).str.len() == 1)
    ]

    after = len(df)
    print(f"Filtered {before - after} indels, {after} SNVs retained", flush=True)

    df.to_csv("mutations_snv.tsv", sep='\\t', index=False)
    EOF
    """
}


// SPLIT MUTATIONS BY PATIENT  (patient mode only)
process SPLIT_MUTATIONS {

    publishDir "${params.outdir}/split_mutations", mode: 'copy'

    input:
    path mutations_tsv

    output:
    path "*.tsv", emit: patients

    script:
    """
    mkdir -p patients

    python3 - << 'EOF'
    import pandas as pd

    df = pd.read_csv("${mutations_tsv}", sep='\\t')

    if 'ID' not in df.columns:
        raise ValueError("mutations TSV must have an 'ID' column for patient mode")

    for pid in df['ID'].dropna().unique():
        safe = str(pid).replace("/", "_").replace(" ", "_")
        sub  = df[df['ID'] == pid].copy()
        sub.to_csv(f"{safe}.tsv", sep='\\t', index=False)
    EOF
    """
}

// TSV -> BED  (normalises chr prefix)
process TO_BED {
    input:
    tuple val(label), path(tsv)

    output:
    tuple val(label), path(tsv), path("${label}.bed")

    script:
    """
    python3 - << 'EOF'
    import pandas as pd

    df = pd.read_csv("${tsv}", sep='\\t')

    df["chrom"] = df["chrom"].astype(str).apply(
        lambda x: x if x.startswith("chr") else "chr" + x
    )

    if df.empty:
        open("${label}.bed", "w").close()
    else:
        bed = pd.DataFrame({
            "chrom": df["chrom"],
            "start": df["position"].astype(int) - 1,
            "end":   df["position"].astype(int)
        }).drop_duplicates()

        bed.to_csv("${label}.bed", sep='\\t', header=False, index=False)
    EOF
    """
}


// Output columns: sample_id, chrom, position, depth, ref_count, alt_count, alt, alt_VAF, origin, bam_file
// Every mutation in the TSV gets a row; positions with no pileup coverage are zero-filled rather than dropped.
process VAF {

    publishDir "${params.outdir}/vaf", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(bai), val(bam_abs_path), path(mutations_tsv), path(bed)
    path reference

    output:
    path "${sample_id}_vaf.tsv"

    script:
    """
    bcftools mpileup \\
        -f ${reference} \\
        -R ${bed} \\
        -B \\
        -d 8000 \\
        -a FORMAT/AD,FORMAT/DP \\
        -Ou ${bam} | \\
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP]\\t[%AD]\\n' \\
        > raw_pileup.tsv

    python3 - << 'EOF'
    import pandas as pd

    mut_df = pd.read_csv("${mutations_tsv}", sep='\\t')
    bam_path = "${bam_abs_path}"

    # make origin optional
    if "origin" not in mut_df.columns:
        if ${params.controls}:
            mut_df["origin"] = "control"
        else:
            mut_df["origin"] = "patient"

    cols = ['chrom','position','ref','alt_list','depth','ad']
    try:
        df = pd.read_csv("raw_pileup.tsv", sep='\\t', names=cols)
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
        chrom = str(mut.chrom)
        pos = int(mut.position)
        origin = str(mut["origin"])
        alt = str(mut["alt"]).upper() if pd.notna(mut["alt"]) else ""
        
        match = bcf_df[(bcf_df.chrom==chrom)&(bcf_df.position==pos)]
        
        if match.empty:
            final.append(["${sample_id}", chrom, pos, 0, 0, 0, alt, 0.0, origin, bam_path])
            continue
        
        row = match.iloc[0]
        alt_count = row.alt_counts[row.alts.index(alt)] if alt in row.alts else 0
        depth = row.depth
        
        final.append(["${sample_id}", chrom, pos, depth, row.ref_count, alt_count, alt, round((alt_count/depth*100) if depth>0 else 0, 4), origin, bam_path])

    out = pd.DataFrame(final, columns=['sample_ID', 'chrom','position','depth','ref_count','alt_count','alt','alt_VAF%', 'origin', 'bam_file'])
    out.to_csv("${sample_id}_vaf.tsv", sep='\\t', index=False)
    EOF
    """
}

workflow {

    mutations_raw = file(params.mutations_tsv)
    reference_ch  = Channel.value(file(params.reference))

    // filter indels first regardless of mode
    snv_tsv = FILTER_MUTATIONS(mutations_raw)

    // CONTROLS MODE  — run every BAM against the full SNV mutations
    if( params.controls ) {

        log.info "RUNNING CONTROLS MODE"

        bam_ch = Channel
            .fromPath("${params.bam_dir}/*.bam")
            .map { bam ->
                def bai = file("${bam}.bai").exists()
                    ? file("${bam}.bai")
                    : file("${bam.parent}/${bam.baseName}.bai")
                tuple(bam.baseName, bam, bai, bam.toAbsolutePath().toString())
            }

        // single BED from the full SNV TSV
        bed_ch = TO_BED( snv_tsv.map { tsv -> tuple("controls", tsv) } )
            .map { label, tsv, bed -> tuple(tsv, bed) }

        // pair every BAM with the shared TSV + BED
        vaf_inputs = bam_ch.combine(bed_ch)
            .map { sample_id, bam, bai, bam_abs, tsv, bed ->
                tuple(sample_id, bam, bai, bam_abs, tsv, bed)
            }

        VAF( vaf_inputs, reference_ch )
    }


    // PATIENT MODE  — split TSV by patient ID, match each to its BAM
    else {

        log.info "RUNNING PATIENT MODE"

        patient_tsvs = SPLIT_MUTATIONS(snv_tsv)
            .flatten()
            .map { tsv ->
                def patient_id = tsv.baseName.replaceFirst("^patient_", "")
                tuple(patient_id, tsv)
            }

        // TO_BED emits: tuple(patient_id, tsv, bed)
        patient_beds = TO_BED(patient_tsvs)

        vaf_inputs = patient_beds
            .map { patient_id, tsv, bed ->
                def bam = file("${params.bam_dir}/${patient_id}.bam")
                def bai = file("${bam}.bai").exists()
                    ? file("${bam}.bai")
                    : file("${bam.parent}/${bam.baseName}.bai")
                tuple(patient_id, bam, bai, bam.toAbsolutePath().toString(), tsv, bed)
            }
            .filter { patient_id, bam, bai, bam_abs, tsv, bed ->
                if( !bam.exists() ) {
                    log.warn "No BAM found for patient ${patient_id} — skipping"
                    return false
                }
                return true
            }

        VAF( vaf_inputs, reference_ch )
    }
}
