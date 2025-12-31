# Coverage assessment; not for ichorCNA, but important
ls *bam | parallel -j 20 --plus \
'mkdir -p /groups/wyattgrp/projects/PANGEN/seq_metrics/mosdepth/{/.}/ && \
mosdepth \
-b "/groups/wyattgrp/projects/PANGEN/hg38.bed" \
-m -t 4 -x -n \
/groups/wyattgrp/projects/PANGEN/seq_metrics/mosdepth/{/.}/{/.} \
{}'

# Generate read counts per 50kb window
ls *bam | parallel -j 20 --plus \
'/groups/wyattgrp/software/hmmcopy_utils/bin/readCounter \
--window 50000 \
--quality 20 \
--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
{} \
> /groups/wyattgrp/projects/PANGEN/ichorCNA/HMMreadCounter/{/.}.wig'
"""
Usage: ./readCounter [options] <BAM file>

Options:
    -w, --window <int>       Specify the size of non-overlapping windows [1000]
    -q, --quality <int>      Specify the mapping quality value below which reads are ignored

    -l, --list               List all chromosomes in BAM reference file
    -s, --sequence <string>  Specify the entries and order of sequences to analyze [ALL],
                             the <string> should be a comma-delimited list (NO spaces)

    -b, --build              Build BAM index for file (same index format as SAMtools)
Example:
    ./readCounter -w 100 -s 1,3,5,X aligned_reads.bam > readcounts.seg
"""


# Make panel of normals from cystoadenoma samples
Rscript /groups/wyattgrp/users/zshong/software/ichorCNA/scripts/createPanelOfNormals.R \
--filelist "/groups/wyattgrp/projects/PANGEN/ichorCNA/normal_wigs.txt" \
--gcWig "/home/jbacon/mambaforge/envs/ichorcna/lib/R/library/ichorCNA/extdata/gc_hg38_50kb.wig" \
--mapWig "/home/jbacon/mambaforge/envs/ichorcna/lib/R/library/ichorCNA/extdata/map_hg38_50kb.wig" \
--centromere "/home/jbacon/mambaforge/envs/ichorcna/lib/R/library/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt" \
--outfile "/groups/wyattgrp/users/zshong/codebook/pangen_PoN.rds"

# Run ichorCNA 
ls *wig | parallel -j 20 --plus \
'Rscript "/groups/wyattgrp/software/ichorCNA/scripts/runIchorCNA.R" \
--id {/.} \
--WIG "/groups/wyattgrp/projects/PANGEN/ichorCNA_standard/HMMreadCounter/{/}" \
--ploidy "c(2,3)" \
--normal "c(0.5,0.6,0.7,0.8,0.9,0.95)" \
--maxCN 5 \
--gcWig "/home/jbacon/mambaforge/envs/ichorcna/lib/R/library/ichorCNA/extdata/gc_hg38_50kb.wig" \
--mapWig "/home/jbacon/mambaforge/envs/ichorcna/lib/R/library/ichorCNA/extdata/map_hg38_50kb.wig" \
--centromere "/home/jbacon/mambaforge/envs/ichorcna/lib/R/library/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt" \
--normalPanel "/groups/wyattgrp/projects/PANGEN/ichorCNA_standard/pangen_PoN.rds_median.rds" \
--includeHOMD False \
--chrs "c(1:22, \"X\")" \
--chrTrain "c(1:22)" \
--estimateNormal True \
--estimatePloidy True \
--estimateScPrevalence True \
--scStates "c(1,3)" \
--txnE 0.9999 \
--txnStrength 10000 \
--outDir "/groups/wyattgrp/projects/PANGEN/ichorCNA_standard/data/"'
