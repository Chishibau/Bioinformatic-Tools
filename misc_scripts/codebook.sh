# wrap command to submit on cluster

sbatch --job-name="cnvkit_${sample_name}" \
           --output="cnvkit_${sample_name}_%j.out" \
           --error="cnvkit_${sample_name}_%j.err" \
           --nobotime=6:00:00 \
           --mem=64G \
           --cpus-per-task=16 \
           --wrap="______"

sbatch --job-name="mutato" \
        --output="/groups/wyattgrp/users/zshong/projects/tmp_outputs/mutato.out" \
        --error="/groups/wyattgrp/users/zshong/projects/tmp_outputs/mutato.err" \
        --time="UNLIMITED" \
        --mem=64G \
        --cpus-per-task=16 \
        --partition=long \
        --wrap="mutato analyze --alt-reads=2 \
               --alt-frac=0.02 \
               --min-mapq=50 \
               --min-end-dist=15 \
               --hotspots="/groups/wyattgrp/users/jbacon/reference/matti_hg38/mutation_hotspots.tsv" \
	       --bg-ratio=0 \
	       --germline-ratio=0 \
               --threads=90 \
               --output-dir="/groups/wyattgrp/projects/uc_metacohort/WGS/ch" \
               "/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa" \
               /groups/wyattgrp/projects/uc_metacohort/WGS/ch/wbc_sample_list.tsv"

sbatch --job-name="maple_copy" \
        --output="/groups/wyattgrp/users/zshong/projects/tmp_outputs/rclone.out" \
        --error="/groups/wyattgrp/users/zshong/projects/tmp_outputs/rclone.err" \
        --time="UNLIMITED" \
        --mem=32G \
        --cpus-per-task=16 \
        --partition=long \
        --wrap="rclone copy --include-from /groups/wyattgrp/users/zshong/projects/hla_wgs/ref/copy_files.txt \
        maple:'cfDNA WGS/NovaSeq X' \
        /groups/wyattgrp/users/zshong/projects/hla_wgs/fqs --progress"
               
python /groups/wyattgrp/users/zshong/codebook/misc_scripts/output_unique_lines.py file1.tsv file2.tsv output.tsv