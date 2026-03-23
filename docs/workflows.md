# Workflows

This page documents all workflows in the pipeline.

## Contents

- [(entry)](#entry) *(entry point)*
- [GRANDEUR](#grandeur)
- [AVERAGE_NUCLEOTIDE_IDENTITY](#average-nucleotide-identity)
- [DE_NOVO_ALIGNMENT](#de-novo-alignment)
- [INITIALIZE](#initialize)
- [PHYLOGENETIC_ANALYSIS](#phylogenetic-analysis)
- [QUALITY_ASSESSMENT](#quality-assessment)
- [REPORT](#report)
- [SUBTYPING](#subtyping)
- [TAXONOMIC_PROFILING](#taxonomic-profiling)
- [TEST](#test)

## (entry) {#entry}

**Entry workflow**

*Defined in `main.nf:26`*


## GRANDEUR {#grandeur}

*Defined in `workflows/grandeur/main.nf:9`*

### Inputs

| Name | Description |
|------|-------------|
| `ch_raw_reads` | - |
| `ch_fastas` | - |
| `ch_reference_genomes` | - |
| `ch_versions` | - |
| `ch_genome_sizes` | - |
| `ch_mash_db` | - |
| `ch_kraken2_db` | - |
| `ch_checkm2_db` | - |
| `ch_sylph_db` | - |
| `dataset_script` | - |
| `evaluat_script` | - |
| `jsoncon_script` | - |
| `multiqc_script` | - |
| `summary_script` | - |
| `summfle_script` | - |
| `version_script` | - |

### Outputs

| Name | Description |
|------|-------------|
| `<none>` | - |


## AVERAGE_NUCLEOTIDE_IDENTITY {#average-nucleotide-identity}

*Defined in `subworkflows/local/average_nucleotide_identity/main.nf:10`*

**Keywords:** `ani`, `nucleotide identity`, `skani`, `reference download`, `taxonomy`, `species estimation`

Runs average nucleotide identity (ANI) analysis. Optionally identifies and downloads  reference genomes from NCBI using datasets, then calculates sketch distances using Skani  to estimate the organism and identify top hits.

### Components

This workflow uses the following modules/subworkflows:

- `spestimator`
- `species`
- `datasets/summary`
- `datasets/download`
- `references`
- `skani/sketch`
- `skani/dist`

### Inputs

| Name | Description |
|------|-------------|
| `meta` | Groovy Map containing sample information e.g. [ id:'test', single_end:false ] |
| `ch_contigs` | Channel containing assembled contigs/fastas to analyze |
| `ch_reference_genomes` | Channel containing user-provided local reference genome fasta files |
| `ch_species` | Channel containing initial species estimations/mash results |
| `dataset_script` | Script required to run the datasets reference genome lookups |

### Outputs

| Name | Description |
|------|-------------|
| `meta` | Groovy Map containing sample information e.g. [ id:'test', single_end:false ] |
| `for_summary` | Combined channel of summary CSV/TSV files from spestimator, datasets, and skani |
| `top_hit` | Tuple containing meta, a list of [species, genus], and the top hit reference fasta file |
| `ch_org_contigs` | Tuple containing meta, predicted organism list [species, genus], and the input contigs |
| `ch_salmonella` | Skani dist output specific to Salmonella hits |
| `ch_ecoli` | Skani dist output specific to E. coli hits |
| `ch_kleb` | Skani dist output specific to Klebsiella hits |
| `ch_gas` | Skani dist output specific to GAS (Group A Strep) hits |
| `ch_strep` | Skani dist output specific to Streptococcus hits |
| `ch_legionella` | Skani dist output specific to Legionella hits |
| `ch_vibrio` | Skani dist output specific to Vibrio hits |
| `ch_acinetobacter` | Skani dist output specific to Acinetobacter hits |
| `ch_myco` | Skani dist output specific to Mycobacterium hits |
| `ch_gc` | Skani dist output specific to Neisseria gonorrhoeae hits |
| `versions` | File containing software versions |

**Authors:** [@erinyoung](https://github.com/erinyoung)


## DE_NOVO_ALIGNMENT {#de-novo-alignment}

*Defined in `subworkflows/local/de_novo_alignment/main.nf:4`*

**Keywords:** `de novo`, `assembly`, `fastp`, `spades`, `fastq`, `alignment`

Runs de novo assembly. Performs FASTQ file filtering using FASTP and de novo assembly using SPAdes.

### Components

This workflow uses the following modules/subworkflows:

- `fastp`
- `spades`

### Inputs

| Name | Description |
|------|-------------|
| `meta` | Groovy Map containing sample information e.g. [ id:'test', single_end:false ] |
| `reads` | Channel containing sample metadata and raw sequencing reads (FASTQ format) |

### Outputs

| Name | Description |
|------|-------------|
| `meta` | Groovy Map containing sample information e.g. [ id:'test', single_end:false ] |
| `reads_contigs` | Tuple containing sample metadata, the cleaned reads, and the assembled contigs |
| `clean_reads` | Channel containing sample metadata and the filtered/cleaned FASTQ files from FASTP |
| `contigs` | Channel containing sample metadata and the assembled contigs from SPAdes |
| `for_multiqc` | FASTP output files (JSON/HTML) to be routed to MultiQC for reporting |
| `versions` | File containing software versions |

**Authors:** [@erinyoung](https://github.com/erinyoung)


## INITIALIZE {#initialize}

*Defined in `subworkflows/local/initialize/main.nf:63`*

**Keywords:** `initialization`, `setup`, `parameters`, `validation`, `configuration`

Initializes the workflow by evaluating parameters, checking inputs, and setting up channels  for reference databases, Python scripts, and input files (FASTQ reads, FASTA assemblies,  and accessions). Optionally runs a test subworkflow to download SRA/genome accessions.

### Components

This workflow uses the following modules/subworkflows:

- `test`

### Inputs

| Name | Description |
|------|-------------|
| `<none>` | - |

### Outputs

| Name | Description |
|------|-------------|
| `reads` | Channel containing sample metadata and paired-end FASTQ files |
| `fastas` | Channel containing sample metadata and input FASTA files |
| `reference_genomes` | Channel containing additional reference genomes for ANI analysis |
| `versions` | File containing software versions (emitted if TEST subworkflow runs) |
| `genome_sizes` | Channel containing the genome sizes JSON file |
| `mash_db` | Channel containing the custom MASH database file |
| `kraken2_db` | Channel containing the Kraken2 database directory |
| `checkm2_db` | Channel containing the CheckM2 database directory |
| `sylph_db` | Channel containing the SYLPH database directory |
| `dataset_script` | Channel containing the datasets_download.py script |
| `evaluat_script` | Channel containing the evaluate.py script |
| `jsoncon_script` | Channel containing the json_convert.py script |
| `multiqc_script` | Channel containing the for_multiqc.py script |
| `summary_script` | Channel containing the summary.py script |
| `summfle_script` | Channel containing the summary_file.py script |
| `version_script` | Channel containing the versions.py script |

**Authors:** [@erinyoung](https://github.com/erinyoung)


## PHYLOGENETIC_ANALYSIS {#phylogenetic-analysis}

*Defined in `subworkflows/local/phylogenetic_analysis/main.nf:13`*

**Keywords:** `phylogeny`, `core genome`, `tree`, `annotation`, `snp distance`, `msa`

Runs phylogenetic analysis. Performs gene prediction (Prokka or Bakta), core genome alignment  (Panaroo or Roary), evaluates core genome robustness, calculates SNP distances (snp-dists),  and constructs/visualizes phylogenetic trees (Mashtree, IQ-TREE, GoTree) and clustered heatmaps.

### Components

This workflow uses the following modules/subworkflows:

- `prokka`
- `bakta`
- `panaroo`
- `roary`
- `core/genome/evaluation`
- `mashtree`
- `ska2`
- `iqtree`
- `gotree`
- `snpdists`
- `heatcluster`

### Inputs

| Name | Description |
|------|-------------|
| `evaluat_script` | Python script (evaluate.py) used for the custom core genome evaluation process |
| `ch_org_contigs` | Tuple containing sample metadata, predicted organism, and assembled contigs |
| `ch_top_hit` | Tuple containing metadata, organism info, and the reference genome from SKANI top hits |

### Outputs

| Name | Description |
|------|-------------|
| `for_multiqc` | Combined channel of output files routed to MultiQC for reporting (e.g., Prokka, Bakta, GoTree, snp-dists, Heatcluster outputs) |
| `summary` | Combined channel of summary TSV/CSV files (e.g., gotree_summary.tsv) |
| `versions` | File containing software versions |

**Authors:** [@erinyoung](https://github.com/erinyoung)


## QUALITY_ASSESSMENT {#quality-assessment}

*Defined in `subworkflows/local/quality_assessment/main.nf:8`*

**Keywords:** `quality control`, `qc`, `fastqc`, `quast`, `amr`, `mlst`, `plasmids`, `checkm2`

Runs quality assessment on raw reads and assembled contigs. Performs read QC with FastQC,  and evaluates assemblies using QUAST, AMRFinderPlus, MLST, PlasmidFinder, and optionally CheckM2.

### Components

This workflow uses the following modules/subworkflows:

- `fastqc`
- `amrfinder`
- `quast`
- `mlst`
- `plasmidfinder`
- `checkm2`

### Inputs

| Name | Description |
|------|-------------|
| `ch_raw_reads` | Channel containing sample metadata and raw paired-end FASTQ reads |
| `ch_clean_reads` | Channel containing sample metadata and cleaned/filtered FASTQ reads |
| `ch_fastas_without_reads` | Channel containing sample metadata and FASTA files that do not have associated reads |
| `ch_all_fastas` | Channel containing sample metadata and all FASTA files (both user-provided and assembled) |
| `ch_reads_contigs` | Tuple containing sample metadata, cleaned reads, and the assembled contigs (used for QUAST) |
| `ch_contigs_org` | Tuple containing sample metadata, organism prediction, and assembled contigs (used for AMRFinderPlus) |
| `ch_checkm2_db` | Channel containing the CheckM2 database directory |
| `summfle_script` | Python script used for summarizing output files |

### Outputs

| Name | Description |
|------|-------------|
| `for_summary` | Combined channel of summary CSV/TSV/JSON files generated by FastQC, AMRFinderPlus, QUAST, MLST, PlasmidFinder, and CheckM2 |
| `for_multiqc` | Combined channel of output files routed to MultiQC for reporting (FastQC and QUAST outputs) |
| `versions` | File containing software versions |

**Authors:** [@erinyoung](https://github.com/erinyoung)


## REPORT {#report}

*Defined in `subworkflows/local/report/main.nf:6`*

**Keywords:** `report`, `multiqc`, `summary`, `versions`, `qc`

Generates final pipeline reports by collating software versions, preparing data for MultiQC,  running MultiQC to create an HTML report, and executing a custom Python script to create  a comprehensive Grandeur summary TSV.

### Components

This workflow uses the following modules/subworkflows:

- `mqc/prep`
- `multiqc`
- `summary`
- `versions`

### Inputs

| Name | Description |
|------|-------------|
| `ch_reads` | Channel containing sample metadata and paired-end FASTQ reads |
| `ch_fastas` | Channel containing sample metadata and FASTA files |
| `for_multiqc` | Combined channel of all output files routed to MultiQC for reporting |
| `for_summary` | Combined channel of summary CSV/TSV/JSON files generated by previous subworkflows |
| `ch_versions` | Channel containing software versions emitted by all executed tools |
| `multiqc_script` | Python script (for_multiqc.py) used to prep outputs before running MultiQC |
| `version_script` | Python script (versions.py) used to convert versions.yml for MultiQC |

### Outputs

| Name | Description |
|------|-------------|
| `summary` | The final generated Grandeur summary TSV file (currently emitted as empty, but tracks the intended SUMMARY output) |
| `versions` | Collated file containing all software versions from the pipeline run |

**Authors:** [@erinyoung](https://github.com/erinyoung)


## SUBTYPING {#subtyping}

*Defined in `subworkflows/local/subtyping/main.nf:17`*

**Keywords:** `subtyping`, `serotyping`, `mlst`, `amr`, `taxonomy`, `typing`

Runs in silico subtyping analysis on genome assemblies. Routes assemblies to specific  tools based on the organism of interest to determine serotypes, sequence types, virulence,  and antimicrobial resistance profiles.

### Components

This workflow uses the following modules/subworkflows:

- `drprg`
- `elgato`
- `emmtyper`
- `json/convert`
- `kaptive`
- `kleborate`
- `meningotype`
- `mykrobe`
- `ngmaster`
- `pbptyper`
- `seqsero2`
- `seqsero2s`
- `serotypefinder`
- `shigapass`

### Inputs

| Name | Description |
|------|-------------|
| `ch_myco` | Channel containing sample metadata and Mycobacterium assemblies |
| `ch_gas` | Channel containing sample metadata and Group A Streptococcus assemblies |
| `ch_kleb` | Channel containing sample metadata and Klebsiella assemblies |
| `ch_legionella` | Channel containing sample metadata and Legionella assemblies |
| `ch_strep` | Channel containing sample metadata and Streptococcus assemblies |
| `ch_salmonella` | Channel containing sample metadata and Salmonella assemblies |
| `ch_ecoli` | Channel containing sample metadata and E. coli / Shigella assemblies |
| `ch_vibrio` | Channel containing sample metadata and Vibrio assemblies |
| `ch_gc` | Channel containing sample metadata and Neisseria gonorrhoeae / meningitidis assemblies |
| `ch_acinetobacter` | Channel containing sample metadata and Acinetobacter assemblies |
| `summfle_script` | Python script used for summarizing output files (used with EMMTYPER and SEROTYPEFINDER) |
| `jsoncon_script` | Python script used for JSON conversion (used with DRPRG) |

### Outputs

| Name | Description |
|------|-------------|
| `for_summary` | Combined channel of summary TSV, CSV, and TXT files generated by the various subtyping tools |
| `versions` | File containing software versions for all executed tools |

**Authors:** [@erinyoung](https://github.com/erinyoung)


## TAXONOMIC_PROFILING {#taxonomic-profiling}

*Defined in `subworkflows/local/taxonomic_profiling/main.nf:6`*

**Keywords:** `taxonomy`, `profiling`, `species identification`, `kraken2`, `mash`, `sylph`, `k-mer`, `minhash`

Rapidly identifies organisms present in sequencing reads and/or assembled FASTA files.  Relies heavily on fast, k-mer and MinHash-based algorithms (Kraken2, Mash, Sylph)  rather than computationally heavy alignments.

### Components

This workflow uses the following modules/subworkflows:

- `kraken2`
- `mash/dist`
- `mash/screen`
- `sylph`

### Inputs

| Name | Description |
|------|-------------|
| `ch_reads` | Channel containing sample metadata and paired-end FASTQ reads |
| `ch_fastas` | Channel containing sample metadata and assembled FASTA files |
| `ch_assemblies` | Channel containing sample metadata and assembled FASTA files (Note/passed into `take` block) |
| `ch_kraken2_db` | Channel containing the Kraken2 database directory |
| `ch_mash_db` | Channel containing the custom MASH database/sketch file |
| `ch_sylph_db` | Channel containing the Sylph database file/directory |

### Outputs

| Name | Description |
|------|-------------|
| `for_ref_download` | Combined channel of species identification results (Kraken2, Mash, Sylph) routed to downstream reference downloading processes |
| `for_summary` | Combined channel of summary CSV/TSV/TXT files generated by Kraken2, Mash, and Sylph |
| `for_multiqc` | Combined channel of output files routed to MultiQC for reporting (Kraken2 outputs) |
| `versions` | File containing software versions for all executed tools |

**Authors:** [@erinyoung](https://github.com/erinyoung)


## TEST {#test}

*Defined in `subworkflows/local/test/main.nf:4`*

**Keywords:** `test`, `download`, `ena`, `ncbi`, `datasets`, `sra`, `fastq`, `fasta`

Downloads test data or specific samples from external databases. Uses ENA's enaDataGet  to download FASTQ files from SRA accessions, and NCBI Datasets to download FASTA  genome assemblies from genome accessions.

### Components

This workflow uses the following modules/subworkflows:

- `ena/download`
- `datasets/download`

### Inputs

| Name | Description |
|------|-------------|
| `ch_sra_accessions` | Channel containing a list of SRA accessions to download from the European Nucleotide Archive (ENA) |
| `ch_genome_accessions` | Channel containing a list of genome accessions to download from NCBI Genomes |

### Outputs

| Name | Description |
|------|-------------|
| `fastq` | Tuple containing sample metadata and paired-end FASTQ reads downloaded from ENA |
| `fasta` | Tuple containing sample metadata and FASTA assembly files downloaded from NCBI |
| `versions` | File containing software versions for the downloading tools |

**Authors:** [@erinyoung](https://github.com/erinyoung)


---

*This pipeline was built with [Nextflow](https://nextflow.io).
Documentation generated by [nf-docs](https://github.com/ewels/nf-docs) v0.2.0 on 2026-03-19 20:49:19 UTC.*
