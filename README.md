# Grandeur Peak

# Warning : What you currently see is under development. It is likely incomplete and incorrect. 

### (We're working on it.)

<img src="https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg" width="250" align="left" />

Named after the beautiful [Grandeur Peak](https://www.alltrails.com/trail/us/utah/grandeur-peak-east-trail-from-church-fork)

Image stolen from [ryancornia](https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg)

Location:  40.707, -111.76, 8,299 ft (2,421 m) summit

More information about the trail leading up to this landmark can be found at [https://utah.com/hiking/grandeur-peak](https://utah.com/hiking/grandeur-peak)

Grandeur Peak is a two-part [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for sequencing of microbial isolates. Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

The "Grandeur" portion of this workflow is intended to be a species agnostic sequencing approach to paired-end Illumina sequencing QC and serotyping in a local public health laboratory. "Grandeur" takes paired-end Illumina reads, cleans them with [seqyclean](https://github.com/ibest/seqyclean), and creates contigs through _de novo_ alignment of the reads with [shovill](https://github.com/tseemann/shovill). Then a variety of QC tools and serotyping tools are employed. This workflow is similar to [Cutshaw](https://staph-b.github.io/staphb_toolkit/workflow_docs/cutshaw/), [Dryad](https://staph-b.github.io/staphb_toolkit/workflow_docs/dryad/), [Foushee](https://staph-b.github.io/staphb_toolkit/workflow_docs/foushee/), and [Tredegar](https://staph-b.github.io/staphb_toolkit/workflow_docs/tredegar/) in the [staphB toolkit](https://github.com/StaPH-B/staphb_toolkit) or [Bactopia](https://github.com/bactopia/bactopia), and the authors or maintainers of this workflow will not be offended an alternative is used.

The "Peak" portion of this workflow uses the resulting contigs in fasta files (or in gff files created by [prokka](https://github.com/tseemann/prokka)) from the "Grandeur" portion through to a phylogenetic tree. The workflow, however, can be used with any fasta - like full genomes downloaded from NCBI, or contigs generated from another workflow. In other words, the end user doesn't need to run "Grandeur" to run "Peak". The gff files created via [prokka](https://github.com/tseemann/prokka), are used by [roary](https://sanger-pathogens.github.io/Roary/) to define and align a core genome (a core genome is the genes that all the provided isolates share, also known as a [pan-genome](https://en.wikipedia.org/wiki/Pan-genome)). This multiple sequence alignment is used to create a tree with [iqtree2](http://www.iqtree.org/) and SNP matrix with [snp-dists](https://github.com/tseemann/snp-dists). As this workflow is dependent on the core genome, all of the files put into this portion of the workflow must be the _same species_ (or predicted to have some sort of similar origin, so related plasmids will work as well). [@erinyoung](https://github.com/erinyoung) strongy recommends that the end user checks 'peak/roary/summary_statistics.txt' to ensure that the number of genes in the core genome makes sense. 

Grandeur Peak will also probably be a workflow of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit) once [@erinyoung](https://github.com/erinyoung) gets around to it and all the containers are ready.

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
- [Docker](https://docs.docker.com/get-docker/) 

# Running Grandeur

```
nextflow run grandeur.nf -c configs/singularity.config
```

"Grandeur" is expecting paired-end Illumina sequences of isolates at `workflow.launchDir + '/reads'`. Keeping in mind that the `workflow.launchDir` is the directory where nextflow command is being run. If the directory of reads is at some other location, it can be specified with `params.reads` in a config file or on the command line as
```
nextflow run grandeur.nf -c configs/singularity.config --reads <directory of reads>
```
A directory will produce files at `workflow.launchDir + '/grandeur'`, but this can also be adjusted with `params.outdir` the same way. 

Additionally, "Graundeur" can use two optional tools to find contamination : [kraken2](https://ccb.jhu.edu/software/kraken2/) or [blobtools](https://blobtools.readme.io/docs) (or both!). These both use large databases that are not included in any container. 

## If the **End User** would like to use [kraken2](https://ccb.jhu.edu/software/kraken2/) to identify contamination:

First, there needs to be a kraken2 database. 

```
mkdir kraken2_db
cd kraken2_db
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz
tar -zxvf minikraken2_v2_8GB_201904.tgz
```
This downloads and expands the directory to 'kraken2_db/minikraken2_v2_8GB_201904_UPDATE'.

- `params.kraken2` needs to be set to `true`
  - `params.kraken2 = 'true'`
- `params.kraken2_db` must be set to the directory with the kraken2 database
  - `params.kraken2_db = '<kraken2 database directory>'`
  - `params.kraken2_db = 'kraken2_db/minikraken2_v2_8GB_201904_UPDATE'` in the example

## If the End User would like to use [blobtools](https://blobtools.readme.io/docs) to identify contamination:

[Blobtools](https://blobtools.readme.io/docs) uses a blast database, so there needs to be a blast database:

```
mkdir blast_db
cd blast_db
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz"
for file in *tar.gz ; do tar -zxvf $file ; done
```
This downloads the NCBI's NT database into 'blast_db'. 
WARNING: The download and decompressing takes a long time. 

- `params.blobtools` needs to be set to `true`
  - `params.blobtools = 'true'`
- `params.blast_db` must be set to the directory with the blast nt database
  - `params.blast_db = '<blast nt database directory'`
  - `params.blast_db = 'blast_db'` in the example

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
│   ├── sample.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt          # Genus and species of the reads
│   └── sample.blobDB.table.txt
├── cg_pipeline
│   ├── sample_cg_pipeline_report.txt                                              # QC metrics of Illumina reads
│   └── cg_pipeline_report.txt
├── contigs
│   └── sample_contigs.fa                                                          # fasta file of contigs
├── fastqc
│   ├── sample_S5_L001_R1_001_fastqc.html
│   ├── sample_S5_L001_R1_001_fastqc.zip
│   ├── sample_S5_L001_R2_001_fastqc.html
│   └── sample_S5_L001_R2_001_fastqc.zip
├── gff
│   └── sample.gff                                                                 # gff file created by prokka
├── grandeur_results.tsv                                                           # summary file
├── kleborate
│   ├── sample_results.txt
│   └── kleborate_results.txt                                                      # klebsiella hypervirulence scoring
├── logs
├── mash
│   ├── sample_mashdist.txt                                                        # mash distances
│   └── sample.msh
├── mlst
│   ├── sample_mlst.txt
│   └── mlst_result.tsv                                                            # mlst of organism (if found)
├── ncbi-AMRFinderplus
├── prokka
│   └── sample
│       ├── sample.err
│       ├── sample.faa
│       ├── sample.ffn
│       ├── sample.fna
│       ├── sample.fsa
│       ├── sample.gbk
│       ├── sample.gff                                                            # annotated contig file that can be used via roary
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
│   └── report.tsv                                                               # QC for contigs
├── seqsero2
│   ├── sample
│   │   ├── sample_H_and_O_and_specific_genes.fasta_mem.fasta
│   │   ├── blasted_output.xml
│   │   ├── data_log.txt
│   │   ├── Extracted_antigen_alleles.fasta
│   │   ├── SeqSero_log.txt
│   │   ├── SeqSero_result.tsv
│   │   └── SeqSero_result.txt
│   └── SeqSero_result.tsv                                                       # Salmonella serotypes
├── seqyclean
│   ├── sample_clean_PE1.fastq.gz                                                # paired-end fastq.gz files after cleaning
│   ├── sample_clean_PE2.fastq.gz
│   ├── sample_clean_SE.fastq.gz
│   ├── sample_clean_SummaryStatistics.tsv
│   ├── sample_clean_SummaryStatistics.txt
│   └── SummaryStatistics.tsv
├── serotypefinder
│   └── sample
│       ├── data.json
│       ├── Hit_in_genome_seq.fsa
│       ├── results_tab.tsv                                                      # E. coli serotypes
│       ├── results.txt
│       ├── Serotype_allele_seq.fsa
│       └── tmp
│           ├── out_H_type.xml
│           └── out_O_type.xml
├── shigatyper
│   └── sample.csv                                                               # Shigatyper serotypes
├── shovill
│   └── sample
│       ├── sample_contigs.fa
│       ├── contigs.fa
│       ├── contigs.gfa
│       ├── shovill.corrections
│       ├── shovill.log
│       └── spades.fasta                                                        # contigs of genome in fasta file
├── shuffled
│   └── sample_shuffled.fastq.gz
└── summary
    ├── sample.summary.tsv
    ├── sample.summary.txt
    └── grandeur_summary.txt                                                    # a table with a summary from all the serotyping and QC tools
```

</details>

# Running Peak

```
nextflow run peak.nf -c configs/singularity.config
```

"Peak" is expecting fasta files, [prokka](https://github.com/tseemann/prokka)'s gff files (which are fasta files that have been annotated with [prokka](https://github.com/tseemann/prokka)), or both. Any file that hasn't already gone through [prokka](https://github.com/tseemann/prokka), will be annotated by [prokka](https://github.com/tseemann/prokka). Then the gff files from [prokka](https://github.com/tseemann/prokka) are aligned with [roary](https://sanger-pathogens.github.io/Roary/). The multiple sequence alignment of the core genome is put through [iqtree2](http://www.iqtree.org/) and [snp-dists](https://github.com/tseemann/snp-dists). 

Fasta files are expected at `workflow.launchDir + '/fastas'`, and can be adjusted with `'params.fastas'`.
GFF files are expected `workflow.launchDir + '/gff'`, and can be adjusted with `'params.gff'`. 

On the command line it looks something like this:
```
nextflow run peak.nf -c configs/singularity.config --fastas <directory of fasta files> --gff <directory of gff files>
```

## Although not required, [@erinyoung](https://github.com/erinyoung) suggests changing some of the default parameters. 
- Changing the prokka parameters to include the genus and species of your organism.  
  - `params.prokka_options` = `'--genus <genus> --species <species>'`
- Changing the iqtree2 parameters to designate which file should be used as an outgroup
  - `params.iqtree2_options` = `'-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000 -o <outgroup>'`

## Additionally, [roary](https://github.com/sanger-pathogens/Roary) can use [kraken2](https://ccb.jhu.edu/software/kraken2/)
[roary](https://github.com/sanger-pathogens/Roary) uses [kraken2](https://ccb.jhu.edu/software/kraken2/) to guess what organisms have been included, and this can be helpful when attempting to determine which files are too different to include. If this feature is desired, the following parameters need to be adjusted :
- `params.kraken2` needs to be set to `true`
  - `params.kraken2 = 'true'`
- `params.kraken2_db` must be set to the directory with the kraken2 database
  - `params.kraken2_db = '<kraken2 database directory>'`

See instructions above for an example way to download a pre-made kraken2 database. 

<details>
   <summary>Final File Tree after running peak.nf</summary>

```
peak/
├── iqtree2
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
│       ├── sample.gff                       # for use in alignment by roary
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
│   └── summary_statistics.txt              # important file with the number of genes involved in core genome
└── snp_dists
    └── snp_matrix.txt                      # SNP matrix counting the number of SNPs that each sample differs by
```

</details>

It'd be nice if a tree was automatically generated from this. Really nice. It has been difficult to find a container that creates a high-quality tree in a command line interface that automatically resizes the text and tree appropriately. (Although, we're always open to suggestions!) The phylogenetic tree found at 'peak/iqtree/iqtree.treefile' or 'peak/iqtree/iqtree.contree' can be visualized through multiple tools, such as [ggtree](https://yulab-smu.top/treedata-book/) or [ete3](http://etetoolkit.org/treeview/).  

# Grandeur Peak wouldn't be possible without:

- [seqyclean](https://github.com/ibest/seqyclean) - cleaning reads
- [shovill](https://github.com/tseemann/shovill) - _de novo_ alignment
- [prokka](https://github.com/tseemann/prokka) - gene annotation - used in peak, optionally in grandeur
- [roary](https://sanger-pathogens.github.io/Roary/) - core genome alignment - only used in peak
- [iqtree2](http://www.iqtree.org/) - phylogenetic tree creation - only used in peak
- [snp-dists](https://github.com/tseemann/snp-dists) - SNP matrix - only used in peak
- [mash](https://github.com/marbl/Mash) - species identifier
- [fastqc](https://github.com/s-andrews/FastQC) - fastq file QC
- [quast](http://quast.sourceforge.net/quast) - contig QC
- [cg-pipeline](https://github.com/lskatz/lyve-SET) - fastq file QC
- [multiqc](https://multiqc.info/) - summarizes QC efforts
- [seqsero2](https://github.com/denglab/SeqSero2) - Salmonella serotyping
- [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper) - Shigella serotyping
- [kleborate](https://github.com/katholt/Kleborate) - Klebsiella serotyping
- [serotypefinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/) - E. coli serotyping
- [mlst](https://github.com/tseemann/mlst) - identification of MLST subtype
- [ncbi-amrfinderplus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) - identification of genes associated with antimicrobial resistence
- [blobtools](https://blobtools.readme.io/docs) - contamination
- [kraken2](https://ccb.jhu.edu/software/kraken2/) - contamination
- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK279684/) - read identification with blobtools
- [bwa](http://bio-bwa.sourceforge.net/) - alignment for blobtools
- [samtools](https://github.com/samtools) - sorting and bam creation for blobtools

# Frequently Asked Questions (aka FAQ)
## What do I do if I encounter an error?

**TELL ME ABOUT IT!!!**
* [Github issue](https://github.com/UPHL-BioNGS/Grandeur/issues)
* [Email me](eriny@utah.gov)
* Send me a message on slack

Be sure to include the command used, what config file was used, and what the **nextflow** error was. 

## Where is an example config file?
There is a template file with all the variables [here](.configs/grandeur_template.config) that the **End User** can copy and edit. All of the parameters are included in that file.

There's also a config file what we use here at UPHL [here](./configs/UPHL.config).

## What about CLIA validation?
At UPHL, we use this workflow to determine the serotype of Salmonella and E. coli under CLIA. Therefore, if you look at [our config file](./configs/UPHL.config), we explicitly specify the containers used for Salmonella serotyping with [seqsero2](https://github.com/denglab/SeqSero2), as well as E. coli serotyping with [serotypefinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/) and [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper). We also specify the containers used prior to these processes in the workflow : [seqyclean](https://github.com/ibest/seqyclean) and [shovill](https://github.com/tseemann/shovill) (which isn't used for [seqsero2](https://github.com/denglab/SeqSero2) or [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper), but it is used for [serotypefinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/)). 

The CLIA officer of the **End User** may request additional locks be put in place, like having all of the containers specified. If additional help is needed, please [submit an issue](https://github.com/UPHL-BioNGS/Grandeur/issues) or [Email me](eriny@utah.gov).

## How were serotyping tools chosen for this workflow?
They perform _well_, their containers were easy to create, and [@erinyoung](https://github.com/erinyoung) had heard about them.

## Are any other tools getting added to "Grandeur"?

As "Grandeur" is intended to be a species agnostic workflow for a local public health laboratory, and sequencing is continuing to expand in its utility, new tools are constantly being needed to analyze isolates to further public health goals. 

For example, the _C. auris_ worfklow [mycosnp](https://github.com/CDCgov/mycosnp) is missing from this workflow, as is any sort of plasmid analysis. 

Many of these additional tools are added by need locally or from the **End User**, so if the **End User** knows of other serotyping/analysis tools, please [submit an issue](https://github.com/UPHL-BioNGS/Grandeur/issues) or tell [@erinyoung](https://github.com/erinyoung) about it, and we'll work in some options. 

[@erinyoung](https://github.com/erinyoung) also appreciates pull requests from forks.

**Warning** : If there's not a relaible container of the suggested tool, [@erinyoung](https://github.com/erinyoung) will request that the **End User** create a container for that tool and contribute to [StaPH-B's docker repositories](https://github.com/StaPH-B/docker-builds).

## What is genome_sizes.json used for?
[genome_sizes.json](./configs/genome_sizes.json) has a list of commonly sequenced organisms and the approximate expected genome size for each organism. This is only used for the "cg-pipeline" process to estimate coverage. A file from the **End User** can be used instead and specified with `params.genome_sizes`. 

## How do I cite this workflow?

This workflow stands on the shoulders of giants. As such, please cite the individual tools that were useful for your manuscript so that those developers can continue to get funding. They are listed above. Mentioning this workflow in the text as "The Grandeur workflow v.<version> (www.github.com/UPHL-BioNGS/Grandeur)" is good enough for [@erinyoung](https://github.com/erinyoung)'s ego.

# Directed Acyclic Diagrams (DAG)
## Full workflow
![alt text](images/grandeur_workflow.png)

![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)


