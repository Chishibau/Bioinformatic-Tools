#!/usr/bin/env python3

import os
import argparse
from pathlib import Path
import textwrap

""" USAGE INFORMATION / README
Generates batch scripts for each FASTQ file pair to convert to SCRAM (lossless FASTQ to CRAM conversion)
NOTE: be aware that different ver of samtools/htslib will output different ver of CRAM file format,
        this script forces all CRAM outputs to be ver3.0

REQUIREMENTS: 
    - have a conda env with samtools and bwa-mem2 installed

    - FASTQ file names need to be in form '{sample}_(1|2).fq*' or {sample}_(1|2).fastq*,
        the sample names must match exactly for the forward and reverse read files

    - have fasta and sam functions installed from https://github.com/annalam/seqkit
        make sure that 'fasta' and 'sam' commands exist within the working environment

    - will create out_dir, slurm array batch script, scram_files, and logs output directories if they don't exist,
        also a samples .tsv file showing what scram files will be generated and corresponding paths
        /abspath/to/out_dir
            ├─ logs
            ├─ scram_files        
            ├─ scram_convert_array.batch        
            ├─ samples.tsv

    - both the reference genome FASTA/index files directory and name (i.e. hg38) needs to be given, 
        ***please ensure that the reference index file is generated using bwa-mem2
        EXAMPLE BELOW:
            /abspath/to/ref_genome_files/
            ├─ hg38.fa
            ├─ hg38.fa.0123         # can have other BWA-MEM2 index files
            ├─ hg38.fa.bwt.2bit.64
            ├─ hg38.fa.pac
            ├─ hg38.fa.ann
            ├─ hg38.fa.amb
            ├─ hg38.fa.sa
            
            In this case the argument values would be:
                ref_name = hg38
                ref_dir = /abspath/to/ref_genome_files/

ARGUMENTS:
    --fq_dir            Directory containing FASTQ files

    --ref_dir           (OPTIONAL) Abs path to directory containing ref genome files
                        default="/groups/wyattgrp/users/jpham/references/hg38/bwamem2_index_matti/"

    --ref_name          (OPTIONAL) Name of reference genome file
                        default="hg38"

    --out_dir           abs path to target dir for script, logs, and converted SCRAM files
                        will automatically make this dir if it doesn't already exist

    --sbatch_time       (OPTIONAL) Time limit for SLURM job, default="00-12:00:00". Provide day-hours:minutes:seconds format.

    --sbatch_cpus       (OPTIONAL) Number of CPU cores for SLURM job, default="8"

    --sbatch_mem        (OPTIONAL) Memory allocation (in GB) for SLURM job, default="32"

    --path_conda        (OPTIONAL) Absolute path to your conda.sh file
                        default=~/anaconda3/etc/profile.d/conda.sh

    --conda_env         (OPTIONAL) Name of existing conda environment with samtools, bwa-mem2 installed
                        default="scram"

Example usage:
    /groups/wyattgrp/users/zshong/projects/to_scram/scram_slurmarray.py \
        --fq_dir /groups/wyattgrp/projects/mannas_collab/fastq \
        --ref_dir /groups/wyattgrp/users/jpham/references/hg38/bwamem2_index_matti/ \
        --ref_name hg38 \
        --out_dir /groups/wyattgrp/users/zshong/projects/to_scram/to_backup \
        --conda_env scram \
        --sbatch_time "UNLIMITED" \
        --sbatch_cpus "16" \
        --sbatch_mem "64"

================================================================      
END USAGE INFORMATION  Author: Zoe Shong
"""

def find_fastq_pairs(fq_dir):
    """ Create nested dictionary for paired-end FASTQ files in the given FASTQ dir
    ================================================================
        Returns:
            samples: {
                "sample_id": {
                    "1": Path("/path/to/sample_1.fastq.gz") or None,
                    "2": Path("/path/to/sample_2.fastq.gz") or None
                },
                ...
            }
    ================================================================
    """
    # find all FASTQ files in given fq_dir
    fq_dir = Path(fq_dir)
    fastqs = sorted(fq_dir.glob("*.fq*")) + sorted(fq_dir.glob("*.fastq*"))

    # build nested dictionary for each sample and corresponding 1 + 2 fq files
    samples = {}
    for fq in fastqs:
        name = fq.name
        if "_1." in name:
            sample_id = name.split("_1.")[0]
            samples.setdefault(sample_id, {}).update({"1": fq})
        elif "_2." in name:
            sample_id = name.split("_2.")[0]
            samples.setdefault(sample_id, {}).update({"2": fq})
        else:
            raise ValueError('Missing or incorrect read direction labeling in file name of:' + name)

    # check if both 1 + 2 fq files exist for every sample
    for sample_id in samples:
        samples[sample_id].setdefault("1", None)
        samples[sample_id].setdefault("2", None)

        if samples[sample_id]["1"] and not samples[sample_id]["2"]:
            raise ValueError('Missing Read 2 FastQ file for: {sample_id}')
        if samples[sample_id]["2"] and not samples[sample_id]["1"]:
            raise ValueError('Missing Read 1 FastQ file for: {sample_id}')

    # check if dictionary of samples is empty
    if not samples:
        raise ValueError('No fq file pairs were found in the given fq_dir')

    return samples


def make_slurm_array_script(
    samples_file, n_samples, dir_ref, ref_name, dir_scram, dir_logs,
    sbatch_time, sbatch_cpus, sbatch_mem, path_conda, conda_env, array_max_concurrent
):

    return textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH --job-name=scram_convert
    #SBATCH --array=1-{n_samples}%{array_max_concurrent}
    #SBATCH --time={sbatch_time}
    #SBATCH --cpus-per-task={sbatch_cpus}
    #SBATCH --mem={sbatch_mem}G
    #SBATCH --output={dir_logs}/scram_%A_%a.log
    #SBATCH --error={dir_logs}/scram_%A_%a.log

    set -euo pipefail

    # Get sample info for this task
    line=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {samples_file})
    sample_id=$(echo $line | cut -f1 -d' ')
    fq1=$(echo $line | cut -f2 -d' ')
    fq2=$(echo $line | cut -f3 -d' ')

    echo "Processing $sample_id"
    echo "FASTQs: $fq1 $fq2"

    source {path_conda}
    conda activate {conda_env}

    scram_file={dir_scram}/$sample_id.scram

    bwa-mem2 mem -t {sbatch_cpus} -v 2 {dir_ref}/{ref_name} $fq1 $fq2 | \\
        samtools view -u -F 2304 | \\
        /groups/wyattgrp/cargo/bin/sam minimize --tags --uncompressed - | \\
        samtools view -C -T {dir_ref}/{ref_name}.fa -O cram,version=3.0 > $scram_file

    fasta_checksum=$(gunzip -c $fq1 $fq2 | /groups/wyattgrp/cargo/bin/fasta checksum --normalize-header -)
    cram_checksum=$(samtools fastq -n $scram_file | /groups/wyattgrp/cargo/bin/fasta checksum --normalize-header -)

    if [[ "$fasta_checksum" != "$cram_checksum" ]]; then
        echo "Checksum mismatch for $sample_id"
        mv $scram_file {dir_scram}/${{sample_id}}_ERROR.scram
        exit 1
    fi

    echo "Completed $sample_id"
    printf "SCRAM conversion successfully completed..."
    """)

def main():
    default_conda_path = os.path.expanduser("~/anaconda3/etc/profile.d/conda.sh")

    parser = argparse.ArgumentParser(
        description="Generate SLURM batch scripts for SCRAM conversion from FASTQ (one script per sample)"
    )

    parser.add_argument("--fq_dir", required=True, help="Directory containing FASTQ files")
    parser.add_argument("--ref_dir", required=False, default="/groups/wyattgrp/users/jpham/references/hg38/bwamem2_index_matti/", help="Absolute path to the directory containing reference genome files and relevant indices.")
    parser.add_argument("--ref_name", required=False, default="hg38", help="Name of reference genome file version.")
    parser.add_argument("--out_dir", required=True, help="Absolute path to directory where all output files will be stored.")
    parser.add_argument("--array_max_concurrent", required=False, default="25", help="Max number of concurrent SLURM array tasks")
    parser.add_argument("--sbatch_time", required=False, default="00-12:00:00", help="Time limit for SLURM job")
    parser.add_argument("--sbatch_cpus", required=False, default="8", help="Number of CPU cores requested for SLURM job")
    parser.add_argument("--sbatch_mem", required=False, default="32", help="Amount of memory (in GB) requested for SLURM job")
    parser.add_argument("--path_conda", required=False, default=default_conda_path, help=f"Absolute path to your conda.sh file")
    parser.add_argument("--conda_env", required=False, default="scram", help=f"Name of existing conda env with tools installed")


    args = parser.parse_args()

    # check if the given output directory exists and create sub directories for output files
    out_dir = Path(args.out_dir)
    dir_scram = out_dir.joinpath('scram_files')
    dir_logs = out_dir.joinpath('logs')

    for d in [out_dir, dir_scram, dir_logs]:
        d.mkdir(parents=True, exist_ok=True)

    # find all fastq paired files in given fq_dir
    samples = find_fastq_pairs(args.fq_dir)

    samples_file = out_dir / "samples.tsv"

    with open(samples_file, "w") as f:
        for sample_id, fq_dict in samples.items():
            f.write(f"{sample_id}\t{fq_dict['1']}\t{fq_dict['2']}\n")


    # generate slurm array batch script for each fastq pair
    array_script = make_slurm_array_script(
                        samples_file=samples_file,
                        n_samples=len(samples),
                        dir_ref=args.ref_dir,
                        ref_name=args.ref_name,
                        dir_scram=dir_scram,
                        dir_logs=dir_logs,
                        sbatch_time=args.sbatch_time,
                        sbatch_cpus=args.sbatch_cpus,
                        sbatch_mem=args.sbatch_mem,
                        path_conda=args.path_conda,
                        conda_env=args.conda_env,
                        array_max_concurrent=args.array_max_concurrent
                        )


    array_script_path = out_dir / "scram_convert_array.batch"
    with open(array_script_path, "w") as f:
        f.write(array_script)

    os.chmod(array_script_path, 0o755)

    print(f"Generated SLURM batch script: {array_script_path}")

    # print output directory absolute paths and instructions for job submission on cluster
    print(f"\nSLURM batch scripts written to: {array_script_path.absolute()}")
    print(f"Logs for batch scripts written to: {dir_logs.absolute()}")
    print(f"SCRAM files will be output to: {dir_scram.absolute()}\n")
    print("Submit script on cluster using:")
    print("============")
    print(f"sbatch {array_script_path.absolute()}")
    print("============\n")


if __name__ == "__main__":
    main()