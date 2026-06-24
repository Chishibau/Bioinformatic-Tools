# region2probe

## Probe Sequence Generation for Tiling of Target Regions

Linux Nextflow pipeline for generating candidate probe sequences from a list of target REGIONS

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
    - numpy
    - biopython
  - bedtools
  - blat
  - seqkit
  - openjdk>=17

**or** alternatively create a new conda environment using the included **conda_env.yml** file:

``` bash
conda env create --name new_env_name --file /path/to/target_dir/probe_from_muts/conda_env.yml
``` 

It is recommended to run this pipeline with a SLURM job executer for faster runtime but can also be ran locally within a Linux environment

### Arguments / Required Files
---

| Argument | Default (if optional) | Description
| --- | ---| --- |
| --regions |  | path to tsv file of mutations
| --ref_fasta |  | reference FASTA file to extract probe sequences
| --probe_size | 120 | length of probe (in bp)
| --blat_threshold | 40 | cutoff threshold for quality of secondary alignment for candidate probe sequences (lower is stricter)
| --bp_homopolymer | 5 | cutoff threshold for homopolymer filtering of probe sequences in bp (lower is stricter, exclusive)
| --output | ./probe_filtering | path to directory for pipeline outputs

#### Required Mutation File Format

Please ensure that the input regions file have the following 6 columns in any order and are tab-separated:
- chr	
- start	
- end	
- id	
- sliding	
- spacing

Only these columns will be used to generate / tile the target regions, any additional columns will NOT be preserved


### Sample Usage / Profiles
---

There are 2 profiles for this pipeline, **slurm** and **local**

Slurm profile is meant to be used on a cluster server with utilizing SLURM job submissions while the local profile should be used in a standard Linux environment

*Slurm profile example:*
``` bash
nextflow /path/to/pipeline/probe_from_muts/probe_filter.nf -profile slurm --ref_fasta /path/to/reference_genome --mutations /path/to/mutations.tsv -bg
```

*Local profile example:*
``` bash
nextflow /path/to/pipeline/probe_from_muts/probe_filter.nf -profile local --ref_fasta /path/to/reference_genome --mutations /path/to/mutations.tsv -bg
```

### Pipeline Output
---

A folder named **region_tiling** will be created in the working directory in which the pipeline was called with the following structure:
```
/path/to/working_dir/
├── .nextflow
├── .nextflow.log.x
├── work/
└── region_tiling/
        ├── filtering/
        |    ├── X_all_ids.txt
        |    ├── X_qual_gc_ids.txt
        |    ├── X_addcrit_invalids.txt
        |    ├── X_homopoly_invalids.txt
        |    └── X_blatI_invalids.txt
        ├── inter_files/
        |    ├── X_cand_probes.bed
        |    ├── X_cand_seqs.fa
        |    └── X_blat_minsI.psl
        ├── qualified/
        |    └── X_qual_probes.txt
        └── results/
            └── X_probes.bed
```

**X** acts as a placeholder for the 'id' value for each row in the input regions.tsv file

#### Output Files Explained
---

**filtering/**

- *X_all_ids.txt* 
: list of all candidate probe_ids for every mutation in the format sample_chr_position::chr:start-end

- *X_qual_gc_ids.txt*
: probe_ids for probe sequences that have a GC content between 30% and 75%

- *X_addcrit_invalids.txt*
: probe_ids with probe sequences that have secondary alignments that either have >3 contiguous perfect matches (≥20 bp each) or blocks >50 bp with >90% bp match

- *X_homopoly_invalids.txt*
: probe_ids with probe sequences that have homopolymers of length greater than threshold defined in arguments

- *X_blatI_invalids*
: probe_ids with probe sequences that have a secondary alignment in the reference genome with a PSL score higher than threshold defined in arguments

**inter_files/**

- *X_cand_probes.bed*
: BED file for all candidate probe locations for all mutations

- *X_cand_seqs.fa*
: FASTA sequences for all candidate probes for all mutations

- *X_blat_minsI.psl*
: BLAT tool output file with PSL score cutoff threshold defined in arguments


**qualified/**

- *X_qual_probes.txt*
: list of all probe_ids that passed all the filtering criterias, may contain more than 1 probe_id for a specific mutation


**results/**

- *X_probes.bed*
: bed file with 4 columns: chr, start, end, probe_id that contains the locations of the generated probes; each file only contains probes for a single given region in the original input .tsv file