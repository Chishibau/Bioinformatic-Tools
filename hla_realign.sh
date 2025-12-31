#!/bin/bash
#SBATCH --job-name=altremap
#SBATCH --output=/groups/wyattgrp/projects/uc_metacohort/WGS/noalt_alignments/logs/altremap_%A_%a.out
#SBATCH --error=/groups/wyattgrp/projects/uc_metacohort/WGS/noalt_alignments/logs/altremap_%A_%a.err
#SBATCH --time=UNLIMITED
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-167%25

# Read the BAM from the list
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /groups/wyattgrp/projects/uc_metacohort/WGS/alignments/bam_paths.txt)
sample=$(basename "$BAM" .bam)

OUTDIR=/groups/wyattgrp/projects/uc_metacohort/WGS/noalt_alignments
REF=/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa
JAR=/groups/wyattgrp/users/zshong/software/bam-tools_v1.4.2.jar

module load java
module load samtools

java -cp $JAR com.hartwig.hmftools.bamtools.remapper.AltContigRemapper \
    -orig_bam_file "$BAM" \
    -ref_genome "$REF" \
    -output_file "${OUTDIR}/${sample}_noalt.bam" \
    -bamtool /groups/share/modules/samtools/1.8/samtools
    -threads 8
