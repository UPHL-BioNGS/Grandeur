# Grandeur Peak

# Warning : What you currently see is under development. It is likely incomplete and incorrect. 

### (We're working on it.)

<img src="https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg" width="250" align="left" />

Named after the beautiful [Grandeur Peak](https://www.alltrails.com/trail/us/utah/grandeur-peak-east-trail-from-church-fork)

Image stolen from [ryancornia](https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg)

Location:  40.707, -111.76, 8,299 ft (2,421 m) summit

More information about the trail leading up to this landmark can be found at [https://utah.com/hiking/grandeur-peak](https://utah.com/hiking/grandeur-peak)

Grandeur Peak is a two-part [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for sequencing of microbial isolates. Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

The "Grandeur" part of the workflow takes paired-end Illumina reads, cleans them with seqyclean, and creates contigs through _de novo_ alignment of the reads with shovill. Then a variety of QC tools and serotyping tools are employed. This workflow, then, is similar to X and X and X in the staphB toolkit, and the authors or maintainers of this workflow will not be offended if you use an alternative.

The "Peak" part of this workflow takes the gff files created via prokka in the "Grandeur" workflow (or any gff files created by prokka), aligns the core genome with roary, and creates a tree with iqtree. This workflow can start with fasta sequences of complete or impartial genomes as well.

Grandeur Peak will also probably be a workflow of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit) once [@erinyoung]("https://github.com/erinyoung") gets around to it and all the containers are ready.

# Getting started

```
git clone https://github.com/UPHL-BioNGS/Grandeur.git
```

To make life easier, follow with

```
cd Grandeur
git init
```

so that `git pull` can be used for updates.

## Prior to starting the workflow

### Install dependencies
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
   - Nextflow version 20+ is required (`nextflow -v` to check)
- [Singularity](https://singularity.lbl.gov/install-linux)

or
- [Docker](https://docs.docker.com/get-docker/) (*with the caveat that the creator and maintainer uses singularity and may not be able to troubleshoot all docker issues*)

# Running `Grandeur`

```
nextflow run grandeur.nf -c configs/singularity.config
```

Grandeur is expecting paired-end Illumina sequences of isolates at `workflow.launchDir + '/reads` and will create a directory with the resulting files at `workflow.launchDir + '/grandeur'`. If the directory of reads is at some other location, it can be specified with `params.reads` in a config file or on the command line as
```
nextflow run grandeur.nf -c configs/singularity.config --reads <directory to reads>
```
The `params.outdir` can also be adjusted accordingly.

<details>
   <summary>Final File Tree after running grandeur.nf</summary>

```
grandeur/
├── aligned
│   ├── sample.sorted.bam
│   └── sample.sorted.bam.csi
├── blastn
│   └── sample.tsv
├── blobtools
│   ├── sample.sample.sorted.bam.cov
│   ├── sample.blobDB.json
│   ├── sample.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png
│   ├── sample.blobDB.json.bestsum.species.p8.span.100.blobplot.read_cov.bam0.png
│   ├── sample.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt
│   └── sample.blobDB.table.txt
├── cg_pipeline
│   ├── sample_cg_pipeline_report.txt
│   └── cg_pipeline_report.txt
├── contigs
│   └── sample_contigs.fa
├── fastqc
│   ├── sample_S5_L001_R1_001_fastqc.html
│   ├── sample_S5_L001_R1_001_fastqc.zip
│   ├── sample_S5_L001_R2_001_fastqc.html
│   └── sample_S5_L001_R2_001_fastqc.zip
├── gff
│   └── sample.gff
├── grandeur_results.tsv
├── kleborate
│   ├── sample_results.txt
│   └── kleborate_results.txt
├── logs
├── mash
│   ├── sample_mashdist.txt
│   └── sample.msh
├── mlst
│   ├── sample_mlst.txt
│   └── mlst_result.tsv
├── ncbi-AMRFinderplus
├── prokka
│   └── sample
│       ├── sample.err
│       ├── sample.faa
│       ├── sample.ffn
│       ├── sample.fna
│       ├── sample.fsa
│       ├── sample.gbk
│       ├── sample.gff
│       ├── sample.log
│       ├── sample.sqn
│       ├── sample.tbl
│       ├── sample.tsv
│       └── sample.txt
├── quast
│   ├── sample
│   │   ├── basic_stats
│   │   │   ├── cumulative_plot.pdf
│   │   │   ├── GC_content_plot.pdf
│   │   │   ├── sample_GC_content_plot.pdf
│   │   │   └── Nx_plot.pdf
│   │   ├── icarus.html
│   │   ├── icarus_viewers
│   │   │   └── contig_size_viewer.html
│   │   ├── quast.log
│   │   ├── report.html
│   │   ├── report.pdf
│   │   ├── report.tex
│   │   ├── report.tsv
│   │   ├── report.txt
│   │   ├── transposed_report.tex
│   │   ├── transposed_report.tsv
│   │   └── transposed_report.txt
│   └── report.tsv
├── seqsero2
│   ├── sample
│   │   ├── sample_H_and_O_and_specific_genes.fasta_mem.fasta
│   │   ├── blasted_output.xml
│   │   ├── data_log.txt
│   │   ├── Extracted_antigen_alleles.fasta
│   │   ├── SeqSero_log.txt
│   │   ├── SeqSero_result.tsv
│   │   └── SeqSero_result.txt
│   └── SeqSero_result.tsv
├── seqyclean
│   ├── sample_clean_PE1.fastq.gz
│   ├── sample_clean_PE2.fastq.gz
│   ├── sample_clean_SE.fastq.gz
│   ├── sample_clean_SummaryStatistics.tsv
│   ├── sample_clean_SummaryStatistics.txt
│   └── SummaryStatistics.tsv
├── serotypefinder
│   └── sample
│       ├── data.json
│       ├── Hit_in_genome_seq.fsa
│       ├── results_tab.tsv
│       ├── results.txt
│       ├── Serotype_allele_seq.fsa
│       └── tmp
│           ├── out_H_type.xml
│           └── out_O_type.xml
├── shigatyper
│   └── sample.csv
├── shovill
│   └── sample
│       ├── sample_contigs.fa
│       ├── contigs.fa
│       ├── contigs.gfa
│       ├── shovill.corrections
│       ├── shovill.log
│       └── spades.fasta
├── shuffled
│   └── sample_shuffled.fastq.gz
└── summary
    ├── sample.summary.tsv
    ├── sample.summary.txt
    └── grandeur_summary.txt
```

</details>

# Running `Peak`

```
nextflow run peak.nf -c configs/singularity.config
```

Grandeur is expecting paired-end Illumina sequences of isolates at `workflow.launchDir + '/reads` and will create a directory with the resulting files at `workflow.launchDir + '/grandeur'`. If the directory of reads is at some other location, it can be specified with `params.reads` in a config file or on the command line as
```
nextflow run peak.nf -c configs/singularity.config --whatever <directory to reads>
```

OR TO CONTIGS


<details>
   <summary>Final File Tree after running grandeur.nf</summary>

```
peak/
├── iqtree
│   ├── iqtree.ckp.gz
│   ├── iqtree.contree                       # treefile without node values
│   ├── iqtree.iqtree
│   ├── iqtree.log
│   ├── iqtree.splits.nex
│   └── iqtree.treefile                      # treefile with node values
├── logs
├── prokka # if starting from fasta
│   └── sample
│       ├── sample.err
│       ├── sample.faa
│       ├── sample.ffn
│       ├── sample.fna
│       ├── sample.fsa
│       ├── sample.gbk
│       ├── sample.gff                       # for use in alignment
│       ├── sample.log
│       ├── sample.sqn
│       ├── sample.tbl
│       ├── sample.tsv
│       └── sample.txt
├── roary
│   ├── accessory_binary_genes.fa
│   ├── accessory_binary_genes.fa.newick
│   ├── accessory_graph.dot
│   ├── accessory.header.embl
│   ├── accessory.tab
│   ├── blast_identity_frequency.Rtab
│   ├── clustered_proteins
│   ├── core_accessory_graph.dot
│   ├── core_accessory.header.embl
│   ├── core_accessory.tab
│   ├── core_alignment_header.embl
│   ├── core_gene_alignment.aln              # core genome alignment
│   ├── fixed_input_files
│   │   └── sample.gff
│   ├── gene_presence_absence.csv
│   ├── gene_presence_absence.Rtab
│   ├── number_of_conserved_genes.Rtab
│   ├── number_of_genes_in_pan_genome.Rtab
│   ├── number_of_new_genes.Rtab
│   ├── number_of_unique_genes.Rtab
│   ├── pan_genome_reference.fa
│   └── summary_statistics.txt.              # important file with the number of genes involved in core genome
└── snp_dists
    └── snp_matrix.txt                       # SNP matrix counting the number of SNPs that each sample differs by
```

</details>

# Grandeur Peak wouldn't be possible without:

- [seqyclean]() - cleaning reads
- [shovill]() - _de novo_ alignment
- [prokka]() - gene annotation
- [roary]() - core genome alignment
- [iqtree]() - phylogenetic tree
- [snp-dists]() - SNP matrix
- [mash]() - species identifier
- [fastqc]() - fastq file QC
- [quast]() - contig QC
- [cg-pipeline]() - fastq file QC
- [seqsero2]() - Salmonella serotyping
- [shigatyper]() - Shigella serotyping
- [kleborate]() - Klebsiella serotyping
- [serotypefinder]() - E. coli serotyping
- [mlst]() - identification of MLST subtype
- [ncbi-amrfinderplus]() - identification of genes associated with antimicrobial resistence
- [blobtools]() - contamination
- [kraken2]() - contamination
- [blastn]() - read identification with blobtools
- [bwa]() - alignment for blobtools
- [samtools]() - sorting and bam creation for blobtools


# Frequently Asked Questions (aka FAQ)
### What do I do if I encounter an error?

**TELL ME ABOUT IT!!!**
* [Github issue](https://github.com/UPHL-BioNGS/Grandeur/issues)
* [Email me](eriny@utah.gov)
* Send me a message on slack

Be sure to include the command used, what config file was used, and what the **nextflow** error was. 

### Where is an example config file?
There is a template file with all the variables [here](.configs/grandeur_template.config) that the **End User** can copy and edit.

There's also a config file what we use here at UPHL [here](./configs/UPHL.config).

### How were serotyping tools chosen for this workflow?
They perform _well_, their containers were easy to create, and []() had heard about them.

If the **End User** knows of other serotyping tools, please let me know and we'll work in some options. 

**Warning** : If there's not a relaible container of the suggested tool, I'll request the **End User** create a container for that tool and contribute to [StaPH-B's docker repositories](https://github.com/StaPH-B/docker-builds).


# Directed Acyclic Diagrams (DAG)
### Full workflow
![alt text](images/grandeur_workflow.png)

![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)


