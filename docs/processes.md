# Processes

This page documents all processes in the pipeline.

## Contents

- [AMRFINDER](#amrfinder)
- [BAKTA](#bakta)
- [CHECKM2](#checkm2)
- [CORE_GENOME_EVALUATION](#core-genome-evaluation)
- [DATASETS_DOWNLOAD](#datasets-download)
- [DATASETS_SUMMARY](#datasets-summary)
- [DRPRG](#drprg)
- [ELGATO](#elgato)
- [EMMTYPER](#emmtyper)
- [ENA_DOWNLOAD](#ena-download)
- [FASTP](#fastp)
- [FASTQC](#fastqc)
- [GOTREE](#gotree)
- [HEATCLUSTER](#heatcluster)
- [IQTREE](#iqtree)
- [JSON_CONVERT](#json-convert)
- [KAPTIVE](#kaptive)
- [KLEBORATE](#kleborate)
- [KRAKEN2](#kraken2)
- [KSNP4](#ksnp4)
- [MASH_DIST](#mash-dist)
- [MASH_SCREEN](#mash-screen)
- [MASHTREE](#mashtree)
- [MENINGOTYPE](#meningotype)
- [MLST](#mlst)
- [MQC_PREP](#mqc-prep)
- [MULTIQC](#multiqc)
- [MYKROBE](#mykrobe)
- [NGMASTER](#ngmaster)
- [PANAROO](#panaroo)
- [PBPTYPER](#pbptyper)
- [PLASMIDFINDER](#plasmidfinder)
- [PROKKA](#prokka)
- [QUAST](#quast)
- [REFERENCES](#references)
- [ROARY](#roary)
- [SEQSERO2](#seqsero2)
- [SEQSERO2S](#seqsero2s)
- [SEROTYPEFINDER](#serotypefinder)
- [SHIGAPASS](#shigapass)
- [SKA2](#ska2)
- [SKANI_DIST](#skani-dist)
- [SKANI_SKETCH](#skani-sketch)
- [SNPDISTS](#snpdists)
- [SOURMASH](#sourmash)
- [SPADES](#spades)
- [SPECIES](#species)
- [SPESTIMATOR](#spestimator)
- [SUMMARY](#summary)
- [SYLPH](#sylph)
- [VERSIONS](#versions)

## AMRFINDER {#amrfinder}

*Defined in `modules/local/amrfinder/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), val(organism), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `amrfinder/*_amrfinder.txt` | `path` | `collect` | - |
| `meta` | `val` | `meta` | - |
| `logs/*/*.log` | `path` | `log` | - |
| `versions.yml` | `path` | `versions` | - |


## BAKTA {#bakta}

*Defined in `modules/local/bakta/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), val(organism), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `bakta/*` | `path` | `bakta_files` | - |
| `bakta/*.txt` | `path` | `for_multiqc` | - |
| `gff/*gff` | `path` | `gff` | - |
| `logs/*/*.log` | `path` | `log` | - |
| `meta` | `val` | `meta` | - |
| `versions.yml` | `path` | `versions` | - |


## CHECKM2 {#checkm2}

*Defined in `modules/local/checkm2/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs), path(db)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("checkm2/*/quality_report.tsv")` | `tuple` | `results` | - |
| `checkm2/*_quality_report.tsv` | `path` | `report` | - |
| `checkm2/*/*` | `path` | `files` | - |
| `logs/${task.process` | `path` | - | - |


## CORE_GENOME_EVALUATION {#core-genome-evaluation}

*Defined in `modules/local/core_genome_evaluation/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `file(fasta), file(summary), file(script)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `core_genome_values.csv` | `path` | `evaluation` | - |
| `core_genome_evaluation/core_genome_evaluation.csv` | `path` | `for_multiqc` | - |
| `logs/${task.process` | `path` | - | - |


## DATASETS_DOWNLOAD {#datasets-download}

*Defined in `modules/local/datasets_download/main.nf:3`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `genomes/*` | `path` | `genomes` | - |
| `versions.yml` | `path` | `versions` | - |


## DATASETS_SUMMARY {#datasets-summary}

*Defined in `modules/local/datasets_summary/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(taxon), file(script)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `datasets/*_genomes.csv` | `path` | `genomes` | - |
| `versions.yml` | `path` | `versions` | - |


## DRPRG {#drprg}

*Defined in `modules/local/drprg/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), val("drprg"), file("drprg/*/*.drprg.json")` | `tuple` | `json` | - |
| `drprg/*/*` | `path` | `results` | - |
| `logs/${task.process` | `path` | - | - |


## ELGATO {#elgato}

*Defined in `modules/local/elgato/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `elgato/*/possible_mlsts.txt` | `path` | `collect` | - |
| `elgato/*/*` | `path` | `results` | - |
| `logs/${task.process` | `path` | - | - |


## EMMTYPER {#emmtyper}

*Defined in `modules/local/emmtyper/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs), file(script)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `emmtyper/*_emmtyper.txt` | `path` | `collect` | - |
| `emmtyper/*` | `path` | `everything` | - |
| `logs/${task.process` | `path` | - | - |


## ENA_DOWNLOAD {#ena-download}

*Defined in `modules/local/ena_download/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(SRR)` | `tuple` | - | - |


## FASTP {#fastp}

*Defined in `modules/local/fastp/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta)` | `tuple` | - | - |


## FASTQC {#fastqc}

*Defined in `modules/local/fastqc/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(fastq)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `fastqc/*html` | `path` | `fastq_files` | - |
| `fastqc/*_fastqc.zip` | `path` | `for_multiqc` | - |
| `fastqc/*_summary.csv` | `path` | `collect` | - |
| `logs/${task.process` | `path` | - | - |


## GOTREE {#gotree}

*Defined in `modules/local/gotree/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `gotree/*.png` | `path` | `for_multiqc` | - |
| `gotree/*` | `path` | `results` | - |
| `gotree/*_stats_all.tsv` | `path` | `stats` | - |
| `logs/${task.process` | `path` | - | - |


## HEATCLUSTER {#heatcluster}

*Defined in `modules/local/heatcluster/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `heatcluster/*` | `path` | `files` | - |
| `heatcluster/*.png` | `path` | `for_multiqc` | - |
| `logs/${task.process` | `path` | - | - |


## IQTREE {#iqtree}

*Defined in `modules/local/iqtree/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `iqtree/iqtree*` | `path` | `tree` | - |
| `iqtree/*nwk` | `path` | `newick` | - |
| `logs/${task.process` | `path` | - | - |


## JSON_CONVERT {#json-convert}

*Defined in `modules/local/json_convert/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), val(analysis), file(json), file(script)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `${analysis` | `path` | - | - |


## KAPTIVE {#kaptive}

*Defined in `modules/local/kaptive/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("kaptive/*.txt")` | `tuple` | `files` | - |
| `kaptive/*_kaptive.tsv` | `path` | `collect` | - |
| `logs/${task.process` | `path` | - | - |


## KLEBORATE {#kleborate}

*Defined in `modules/local/kleborate/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contig), file(script)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `kleborate/*_kleborate.tsv` | `path` | `collect` | - |
| `val(meta), file("kleborate/*_output.txt")` | `tuple` | `result` | - |
| `logs/${task.process` | `path` | - | - |


## KRAKEN2 {#kraken2}

*Defined in `modules/local/kraken2/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(fastq), path(kraken2_db)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `kraken2/*_kraken2_report.txt` | `path` | `for_multiqc` | - |
| `kraken2/*` | `path` | `files` | - |
| `logs/${task.process` | `path` | - | - |


## KSNP4 {#ksnp4}

*Defined in `modules/local/ksnp4/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `ksnp4/*` | `path` | `files` | - |
| `ksnp4/tree.parsimony.tre` | `path` | `newick` | - |
| `ksnp4/tree.ML.tre` | `path` | `tree_ml` | - |
| `ksnp4/SNPs_all_matrix.fasta` | `path` | `snp_matrix` | - |
| `logs/${task.process` | `path` | - | - |


## MASH_DIST {#mash-dist}

*Defined in `modules/local/mashdist/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads), file(reference)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("mash/*.mashdist.txt")` | `tuple` | `mashdist` | - |
| `mash/*_mashdist_summary.csv` | `path` | `results` | - |
| `*err` | `path` | `mash_err` | - |
| `logs/${task.process` | `path` | - | - |


## MASH_SCREEN {#mash-screen}

*Defined in `modules/local/mashscreen/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads), file(reference)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("mash/*.mashscreen.txt")` | `tuple` | `screen` | - |
| `mash/*_mashscreen_summary.csv` | `path` | `screen_results` | - |
| `versions.yml` | `path` | `versions` | - |


## MASHTREE {#mashtree}

*Defined in `modules/local/mashtree/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `mashtree/*` | `path` | `tree` | - |
| `mashtree/*.nwk` | `path` | `newick` | - |
| `logs/${task.process` | `path` | - | - |


## MENINGOTYPE {#meningotype}

*Defined in `modules/local/meningotype/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("meningotype/*.tsv")` | `tuple` | `files` | - |
| `*meningotype.tsv` | `path` | `summary` | - |
| `versions.yml` | `path` | `versions` | - |
| `meta` | `val` | `meta` | - |


## MLST {#mlst}

*Defined in `modules/local/mlst/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contig)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("mlst/*_mlst.txt")` | `tuple` | `files` | - |
| `*_mlst_summary.txt` | `path` | `collect` | - |
| `versions.yml` | `path` | `versions` | - |


## MQC_PREP {#mqc-prep}

*Defined in `modules/local/mqc_prep/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `*mqc*` | `path` | `for_multiqc` | - |


## MULTIQC {#multiqc}

*Defined in `modules/local/multiqc/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `multiqc/multiqc_report.html` | `path` | `report` | - |
| `multiqc/multiqc_data/*` | `path` | `data_folder` | - |
| `logs/${task.process` | `path` | - | - |


## MYKROBE {#mykrobe}

*Defined in `modules/local/mykrobe/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `mykrobe/*.csv` | `path` | `collect` | - |
| `mykrobe/*.json` | `path` | `json` | - |
| `logs/${task.process` | `path` | - | - |


## NGMASTER {#ngmaster}

*Defined in `modules/local/ngmaster/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("ngmaster/*")` | `tuple` | `files` | - |
| `*ngmaster.csv` | `path` | `collect` | - |
| `logs/${task.process` | `path` | - | - |


## PANAROO {#panaroo}

*Defined in `modules/local/panaroo/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `panaroo/*` | `path` | `files` | - |
| `path("panaroo/core_gene_alignment.aln"), path("panaroo/gene_presence_absence.Rtab")` | `tuple` | `core_gene_alignment` | - |
| `logs/${task.process` | `path` | - | - |


## PBPTYPER {#pbptyper}

*Defined in `modules/local/pbptyper/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `pbptyper/${meta.id` | `path` | - | - |


## PLASMIDFINDER {#plasmidfinder}

*Defined in `modules/local/plasmidfinder/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(file)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("plasmidfinder/*/*")` | `tuple` | `files` | - |
| `plasmidfinder/*/*json` | `path` | `collect` | - |
| `logs/${task.process` | `path` | - | - |


## PROKKA {#prokka}

*Defined in `modules/local/prokka/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), val(organism), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `prokka/*/*` | `path` | `prokka_files` | - |
| `prokka/*/*.txt` | `path` | `for_multiqc` | - |
| `gff/*.gff` | `path` | `gff` | - |
| `logs/${task.process` | `path` | - | - |


## QUAST {#quast}

*Defined in `modules/local/quast/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs), file(reads)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `quast/*` | `path` | `files` | - |
| `quast/*_quast_report.tsv` | `path` | `for_multiqc` | - |
| `val(meta), file("quast/*_quast_report.tsv")` | `tuple` | `results` | - |
| `quast/*quast_transposed_report.tsv` | `path` | `collect` | - |
| `quast/*quast_transposed_report_contig.tsv` | `path` | `collect_contig` | - |
| `logs/${task.process` | `path` | - | - |


## REFERENCES {#references}

*Defined in `modules/local/references/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `ref/*` | `path` | `fastas` | - |


## ROARY {#roary}

*Defined in `modules/local/roary/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `roary/*` | `path` | `files` | - |
| `path("roary/core_gene_alignment.aln"), path("roary/gene_presence_absence.Rtab")` | `tuple` | `core_gene_alignment` | - |
| `logs/${task.process` | `path` | - | - |


## SEQSERO2 {#seqsero2}

*Defined in `modules/local/seqsero2/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(file)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("seqsero2/*/*")` | `tuple` | `files` | - |
| `seqsero2/*_seqsero_result.tsv` | `path` | `collect` | - |
| `logs/${task.process` | `path` | - | - |


## SEQSERO2S {#seqsero2s}

*Defined in `modules/local/seqsero2s/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(file)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("seqsero2s/*/*")` | `tuple` | `files` | - |
| `seqsero2s/*_seqsero2s_result.tsv` | `path` | `collect` | - |
| `logs/${task.process` | `path` | - | - |


## SEROTYPEFINDER {#serotypefinder}

*Defined in `modules/local/serotypefinder/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(file), file(script)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `serotypefinder/*/*` | `path` | `files` | - |
| `serotypefinder/*_serotypefinder.tsv` | `path` | `collect` | - |
| `logs/${task.process` | `path` | - | - |


## SHIGAPASS {#shigapass}

*Defined in `modules/local/shigapass/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `shigapass/*_shigapass.tsv` | `path` | `summary` | - |
| `val(meta), file("shigapass/*/*")` | `tuple` | `all_files` | - |
| `logs/${task.process` | `path` | - | - |


## SKA2 {#ska2}

*Defined in `modules/local/ska2/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `ska/*_alignment.aln` | `path` | `aln` | - |
| `ska/*` | `path` | `files` | - |
| `logs/${task.process` | `path` | - | - |


## SKANI_DIST {#skani-dist}

*Defined in `modules/local/skanidist/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file(contigs), file("skani/*.tsv")` | `tuple` | `hits` | - |
| `val(meta), file("skani/*.tsv")` | `tuple` | `results` | - |
| `val(meta), file("Mycobacteri/*")` | `tuple` | `myco` | - |
| `val(meta), file("gas/*")` | `tuple` | `gas` | - |
| `val(meta)` | `tuple` | - | - |


## SKANI_SKETCH {#skani-sketch}

*Defined in `modules/local/skanisketch/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `skani_db` | `path` | `db` | - |
| `logs/${task.process` | `path` | - | - |


## SNPDISTS {#snpdists}

*Defined in `modules/local/snpdists/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `snp-dists/snpdists*.txt` | `path` | `snp_matrix` | - |
| `versions.yml` | `path` | `versions` | - |


## SOURMASH {#sourmash}

*Defined in `modules/local/sourmash/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads), file(reference)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `sourmash/*.search.csv` | `path` | `search` | - |
| `val(meta), file("sourmash/*.summary.csv")` | `tuple` | `results` | - |
| `logs/${task.process` | `path` | - | - |


## SPADES {#spades}

*Defined in `modules/local/spades/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `spades/*/*` | `path` | `files` | - |
| `val(meta), file("contigs/*_contigs.fa")` | `tuple` | `contigs` | - |
| `val(meta), file("contigs/*_contigs.fa"), file(reads)` | `tuple` | `reads_contigs` | - |
| `logs/${task.process` | `path` | - | - |


## SPECIES {#species}

*Defined in `modules/local/species/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `datasets/species_list.txt` | `path` | `species` | - |
| `datasets/accessions_list.txt` | `path` | `accessions` | - |


## SPESTIMATOR {#spestimator}

*Defined in `modules/local/spestimator/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(contigs)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("spestimator/*.csv")` | `tuple` | `results` | - |
| `logs/${task.process` | `path` | - | - |


## SUMMARY {#summary}

*Defined in `modules/local/summary/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `grandeur_summary.tsv` | `path` | `summary_tsv` | - |
| `grandeur_summary.txt` | `path` | `summary_txt` | - |
| `summary/grandeur_extended_summary.tsv` | `path` | `extended_tsv` | - |
| `summary/grandeur_extended_summary.txt` | `path` | `extended_txt` | - |


## SYLPH {#sylph}

*Defined in `modules/local/sylph/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads), path(db)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("sylph/*sylph.tsv")` | `tuple` | `tsv` | - |
| `sylph/*sylph_results.tsv` | `path` | `results` | - |
| `download/*.txt` | `path` | `for_download` | - |
| `logs/${task.process` | `path` | - | - |


## VERSIONS {#versions}

*Defined in `modules/local/versions/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `software_versions_mqc.yml` | `path` | `for_multiqc` | - |
| `software_versions.yml` | `path` | `yml` | - |


---

*This pipeline was built with [Nextflow](https://nextflow.io).
Documentation generated by [nf-docs](https://github.com/ewels/nf-docs) v0.2.0 on 2026-03-19 20:49:19 UTC.*
