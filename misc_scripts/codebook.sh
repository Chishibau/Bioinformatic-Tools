# wrap command to submit on cluster

sbatch --job-name="cnvkit_${sample_name}" \
           --output="cnvkit_${sample_name}_%j.out" \
           --error="cnvkit_${sample_name}_%j.err" \
           --time=6:00:00 \
           --mem=64G \
           --cpus-per-task=16 \
           --wrap="______"

sbatch --job-name="cnvkit_batch" \
        --output="/groups/wyattgrp/projects/brachy_tracks/cnvkit/min25kb/cnvkit_batch.out" \
        --error="/groups/wyattgrp/projects/brachy_tracks/cnvkit/min25kb/cnvkit_batch.err" \
        --time="UNLIMITED" \
        --mem=64G \
        --cpus-per-task=16 \
        --wrap="cnvkit.py batch -m wgs -r /groups/wyattgrp/projects/brachy_tracks/cnvkit/ref/flat_reference_min25kb.cnn -p 10 -d /groups/wyattgrp/projects/brachy_tracks/cnvkit/min25kb /groups/wyattgrp/projects/brachy_tracks/cram_alignments/*FFPE*.cram"


python /groups/wyattgrp/users/zshong/codebook/misc_scripts/output_unique_lines.py file1.tsv file2.tsv output.tsv


