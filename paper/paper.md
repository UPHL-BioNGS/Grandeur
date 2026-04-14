---
title: 'Grandeur: A modular, species-agnostic Nextflow pipeline for microbial genome assembly, typing, and phylogenetics'
tags:
  - nextflow
  - bioinformatics
  - genomics
  - microbiology
  - public health
  - whole genome sequencing
authors:
  - name: Erin Young
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: Utah Public Health Laboratory (UPHL)
   index: 1
date: 14 April 2026
bibliography: paper.bib
---

# Summary

**Grandeur** is a species-agnostic sequencing analysis workflow developed at the Utah Public Health Laboratory (UPHL). Built on Nextflow, the pipeline provides quality control, *de novo* assembly, taxonomic profiling, and in silico serotyping for paired-end Illumina data. 

Grandeur is modular and executes several distinct stages. The pipeline performs *de novo* alignment to clean reads with fastp and assemble them with SPAdes. It features rapid taxonomic profiling via SKANI, Kraken2, Mash, and Sylph, as well as assembly quality assessment using QUAST, CheckM2, and FastQC. Additionally, it supports organism-specific subtyping using tools like Kleborate, SeqSero2, and Legionella SBT, and optionally conducts phylogenetic analysis through core genome alignment and maximum likelihood tree generation. 

# Statement of need

Public health laboratories frequently process diverse, novel, or unknown microbial isolates that require a highly flexible, species-agnostic analytical approach. While intended to augment the CDC's PHOENIX workflow, Grandeur also functions as a powerful standalone pipeline. It is explicitly designed as a short-read *de novo* assembly pipeline with serotyping capabilities. 

By unifying these disparate tools—from read cleaning to phylogenetic tree creation—into a single Nextflow architecture, Grandeur reduces the bioinformatics bottleneck for microbiologists. It enables reproducible execution via containerization and standardizes output reporting, allowing researchers to quickly transform raw FASTQ reads or existing FASTA assemblies into actionable genomic and epidemiological insights.

# Acknowledgements

Grandeur was developed by Erin Young at the Utah Public Health Laboratory (UPHL). Grandeur wouldn't be possible without the open-source bioinformatics tools it wraps, including nf-core tools for schema functionality and nf-docs for documentation generation. We also acknowledge the developers of the foundational software utilized throughout the pipeline.

# References
