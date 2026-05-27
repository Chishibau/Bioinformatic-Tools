## Probe Sequence Generation from Mutations

Linux Nextflow pipeline for calculating VAF (variant allele frequency) for BAM files from a list of target mutations

### Installation
---

Clone repository at desired path in Linux environment with:

``` bash
wget -P /path/to/target_dir https://github.com/Chishibau/Bioinformatic-Tools.git
```

### Dependencies / Conda environment
---

Ensure the following conda packages are installed in the working environment:
  - nextflow
  - python>=3.10
    - pandas
  - bcftools

**or** alternatively create a new conda environment using the included **conda_env.yml** file:

``` bash
conda env create --name new_env_name --file /path/to/target_dir/mrd_seq_processing/conda_env.yml
``` 

It is recommended to run this pipeline with a SLURM job executer for faster runtime but can also be ran locally within a Linux environment

### Arguments / Required Files
---

| Argument | Default (if optional) | Description
| --- | ---| --- |
| --mutations_tsv |  | path to tsv file of mutations
| --bam_dir |  | path to directory containing mutation tsv files
| --reference | "/groups/wyattgrp/reference/hg38/hg38.fa" | reference FASTA file to extract probe sequences
| --controls | false | number of basepairs to consider on each side of mutation
| --outdir | ./results | length of probe (in bp)

#### Required Mutation File Format

Please ensure that the input mutations tsv file has the following columns in any order and are tab-separated:
- ID (optional for controls)
- chrom
- position
- ref
- alt
- origin (optional, defaults to controls/patient depending on input arguments)

Example mutations TSV:
```
chrom   ID	        position  origin	ref	alt
chr8    BT-015-FFPE	4886993   patient	T	C
chr5	BT-011-FFPE	5447675	  patient	G	T
```

### Control vs. Non-control mode

The --controls parameter determines if the pipeline will run VAF calculations for all mutations against ALL BAM files in the given directory or attempt to match sample-specific mutations to the corresponding BAM file.

Control mode is good for variants meant to be used for determining limit of detection and sequencing metrics.

Non-control mode is essentially patient-specific where the ID column of the mutations TSV will be used to find a matching BAM file.

### Sample Usage / Profiles
---

There are 2 profiles for this pipeline, **slurm** and **local**

Slurm profile is meant to be used on a cluster server with utilizing SLURM job submissions while the local profile should be used in a standard Linux environment

*Slurm profile example:*
``` bash
nextflow /path/to/pipeline/mrd_processing.nf --mutations_tsv /path/to/mutations.tsv --bam_dir /path/to/bam_alignments -profile slurm --controls true -bg
```

*Local profile example:*
``` bash
nextflow /path/to/pipeline/mrd_processing.nf --mutations_tsv /path/to/mutations.tsv --bam_dir /path/to/bam_alignments -profile local --controls false -bg
```

### Pipeline Output
---

A folder named **vaf** will be created in the working directory in which the pipeline was called which contains .tsv files of the calcuated VAFs for the processed BAMs:

Sample output .tsv for *non-controls* mode:
```
sample_ID	chrom	position	depth	ref_count	alt_count	alt	alt_VAF%	origin	bam_file
BT-004-FFPE	chr20	38504188	69	55	11	C	15.942	patient	/groups/wyattgrp/projects/brachy_tracks/bam_alignments/BT-004-FFPE.bam
BT-004-FFPE	chr15	35340966	78	61	17	C	21.7949	patient	/groups/wyattgrp/projects/brachy_tracks/bam_alignments/BT-004-FFPE.bam
```

Sample output .tsv for *controls* mode:
```
sample_ID	chrom	position	depth	ref_count	alt_count	alt	alt_VAF%	origin	bam_file
CTRL-VAF-01_BT-MRD	chr4	1794007	1269	1266	3	A	0.2364	TWIST_std	/groups/wyattgrp/projects/mrd_panel_test/alignments/CTRL-VAF-01_BT-MRD.bam
CTRL-VAF-01_BT-MRD	chr4	1799784	3503	3498	5	T	0.1427	TWIST_std	/groups/wyattgrp/projects/mrd_panel_test/alignments/CTRL-VAF-01_BT-MRD.bam
```