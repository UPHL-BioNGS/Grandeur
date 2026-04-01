# Pipeline Inputs

This page documents all input parameters for the pipeline.

## Input/output options

### `--input` {#input}

**Type:** `string` | *Optional* | **Format:** `file-path`

Path to comma-separated file containing information about the samples in the experiment.

> You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.

**Pattern:** `^\S+\.csv$`


### `--sample_sheet` {#sample-sheet}

**Type:** `string` | *Optional*

csv with sample,read1,read2


### `--fastas` {#fastas}

**Type:** `string` | *Optional*

directory with fasta files (not compatible with cloud resources)


### `--fasta_list` {#fasta-list}

**Type:** `string` | *Optional*

A sample sheet for fasta files


### `--reads` {#reads}

**Type:** `string` | *Optional*

directory with paired-end illumina fastq files (not compatible with cloud resources)


### `--outdir` {#outdir}

**Type:** `string` | **Required** | **Format:** `directory-path`

The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.

**Default:** `grandeur`


### `--sra_accessions` {#sra-accessions}

**Type:** `string` | *Optional*

list of accessions to download from the SRA

**Default:** `[]`


### `--genome_accessions` {#genome-accessions}

**Type:** `string` | *Optional*

list of accessions to download from genomes

**Default:** `[]`


### `--genome_sizes` {#genome-sizes}

**Type:** `string` | *Optional*

text of genome sizes



## Reference files/paths

### `--checkm2_db` {#checkm2-db}

**Type:** `string` | *Optional*

prepared checkm2 reference file


### `--kraken2_db` {#kraken2-db}

**Type:** `string` | *Optional*

directory of kraken2 database


### `--mash_db` {#mash-db}

**Type:** `string` | *Optional*

prepared mash reference msh file


### `--sylph_db` {#sylph-db}

**Type:** `string` | *Optional*

prepared sylph reference file


### `--reference_genomes` {#reference-genomes}

**Type:** `string` | *Optional*

list of genomes (in fasta format) for ANI references



## workflow values

### `--min_core_genes` {#min-core-genes}

**Type:** `integer` | *Optional*

minimum number of genes in core genome alignment for iqtree2 (default is 500)

**Default:** `500`


### `--min_core_per` {#min-core-per}

**Type:** `number` | *Optional*

minimum percentage number of core genes in core genome alignment for iqtree2 (default is 0.5 or 50%)

**Default:** `0.5`


### `--minimum_reads` {#minimum-reads}

**Type:** `integer` | *Optional*

the minimum number of reads in a fastq file required to move to de novo alignment

**Default:** `10000`



## Subworkflow toggles

### `--exclude_top_hit` {#exclude-top-hit}

**Type:** `boolean` | *Optional*

removes ANI top hit from msa


### `--msa` {#msa}

**Type:** `boolean` | *Optional*

toggles whether or not phylogenetic analysis will be run on samples


### `--aligner` {#aligner}

**Type:** `string` | *Optional*

chooses core genome aligner (params.msa must be set to true)

**Default:** `panaroo`

**Allowed values:**
- `roary`
- `panaroo`


### `--annotator` {#annotator}

**Type:** `string` | *Optional*

chooses annotator (params.msa must be set to true)

**Default:** `bakta`

**Allowed values:**
- `bakta`
- `prokka`


### `--skip_extras` {#skip-extras}

**Type:** `boolean` | *Optional*

turns off blobtools, kraken2, skani, mash, and report generation subworkflows


### `--current_datasets` {#current-datasets}

**Type:** `boolean` | *Optional*

toggles whether or not genomes are downloaded from NCBI



## Generic options

### `--help` {#help}

**Type:** `boolean` | *Optional*

Display help text.


### `--version` {#version}

**Type:** `boolean` | *Optional*

Display version and exit.


### `--publish_dir_mode` {#publish-dir-mode}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.

**Default:** `copy`


### `--email_on_fail` {#email-on-fail}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.


### `--plaintext_email` {#plaintext-email}

**Type:** `boolean` | *Optional*

Stolen from example and might not do anything.


### `--monochrome_logs` {#monochrome-logs}

**Type:** `boolean` | *Optional*

Stolen from example and might not do anything.


### `--email` {#email}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.


### `--hook_url` {#hook-url}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.


### `--pipelines_testdata_base_path` {#pipelines-testdata-base-path}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.

**Default:** `https://raw.githubusercontent.com/nf-core/test-datasets/`


### `--config_profile_description` {#config-profile-description}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.


### `--config_profile_name` {#config-profile-name}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.


### `--custom_config_version` {#custom-config-version}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.

**Default:** `master`


### `--custom_config_base` {#custom-config-base}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.

**Default:** `https://raw.githubusercontent.com/nf-core/configs/master`


### `--config_profile_url` {#config-profile-url}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.


### `--config_profile_contact` {#config-profile-contact}

**Type:** `string` | *Optional*

Stolen from example and might not do anything.



---

*This pipeline was built with [Nextflow](https://nextflow.io).
Documentation generated by [nf-docs](https://github.com/ewels/nf-docs) v0.2.0 on 2026-03-19 20:49:19 UTC.*
