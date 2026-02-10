## Probe Sequence Generation from Mutations

Linux Nextflow pipeline for generating candidate probe sequences from a list of target mutations

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
| --mutations |  | path to tsv file of mutations
| --mutations_dir |  | path to directory containing mutation tsv files
| --ref_fasta |  | reference FASTA file to extract probe sequences
| --size | 100 | number of basepairs to consider on each side of mutation
| --probe_size | 120 | length of probe (in bp)
| --blat_threshold | 40 | cutoff threshold for quality of secondary alignment for candidate probe sequences (lower is stricter)
| --bp_homopolymer | 5 | cutoff threshold for homopolymer filtering of probe sequences in bp (lower is stricter, exclusive)

#### Required Mutation File Format

Please ensure that the input mutation files have the following 3 columns in any order and are tab-separated:
- Sample
- Chrom
- Position

Only these columns will be used to generate unique probe_ids for each mutation, but the input files may contain additional columns that will be preserved in the final output file


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

A folder named **probe_filtering** will be created in the working directory in which the pipeline was called with the following structure:
```
/path/to/working_dir/
├── .nextflow
├── .nextflow.log.x
├── work/
└── probe_filtering/
    └── mutations_filename/
        ├── filtering/
        |    ├── all_probe_ids.txt
        |    ├── filtered_gc_ids.txt
        |    ├── addcriteria_invalid_ids.txt
        |    ├── homopolymers_invalid_ids.txt
        |    └── secondary40_invalid-probes.psl
        ├── inter_files/
        |    ├── all_cand_probes.bed
        |    ├── all_cand_seqs.fa
        |    ├── all_muts_Xbp_regions.bed
        |    └── blat_all_minsX.psl
        ├── qualified/
        |    ├── probe_loc_stats.tsv
        |    └── qualified_probes.txt
        └── results/
            ├── muts_with_probe_label.tsv
            ├── close_loc_muts.tsv
            └── no_close_muts.tsv
```

If pipineline is called with multiple mutation files, each mutation file will receive its own output directory within **probe_filtering**

#### Output Files Explained
---

**filtering/**

- *all_probe_ids.txt* 
: list of all candidate probe_ids for every mutation in the format sample_chr_position::chr:start-end

- *filtered_gc_ids.txt*
: probe_ids for probe sequences that have a GC content between 30% and 75%

- *addcriteria_invalid_ids.txt*
: probe_ids with probe sequences that have secondary alignments that either have >3 contiguous perfect matches (≥20 bp each) or blocks >50 bp with >90% bp match

- *homopolymers_invalid_ids.txt*
: probe_ids with probe sequences that have homopolymers of length greater than threshold defined in arguments

- *secondaryX_invalid-probes.psl*
: probe_ids with probe sequences that have a secondary alignment in the reference genome with a PSL score higher than threshold defined in arguments

**inter_files/**

- *all_cand_probes.bed*
: BED file for all candidate probe locations for all mutations

- *all_cand_seqs.fa*
: FASTA sequences for all candidate probes for all mutations

- *all_muts_Xbp_regions.bed*
: BED file with window sizes defined in arguments centered around each mutation

- *blat_all_minsX.psl*
: BLAT tool output file with PSL score cutoff threshold defined in arguments


**qualified/**

- *probe_loc_stats.tsv*
: table of qualified probe_ids with mutation closest to center with columns indicating distance of mutation from the start of the probe sequence

- *qualified_probes.txt*
: list of all probe_ids that passed all the filtering criterias, may contain mroe than 1 probe_id for a specific mutation


**results/**

- *muts_with_probe_label.tsv*
: original input mutation file with added column of probe_id, column value will be empty if there are no valid probes for a given mutation

- *close_loc_muts.tsv*
: table of mutations that have a distance of less than half the length of a probe between them, can be used to determined if 1 probe can be used for more than 1 mutation

- *no_close_muts.tsv*
: produced by removing the lines of close_loc_muts.tsv from the muts_with_probe_label.tsv file