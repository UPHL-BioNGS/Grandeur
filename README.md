# Grandeur Peak

<img src="https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg" width="500" align="left" />

Named after the beautiful [Grandeur Peak](https://www.alltrails.com/trail/us/utah/grandeur-peak-east-trail-from-church-fork)

Image Credit: [ryancornia](https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg)

Location:  40.707, -111.76, 8,299 ft (2,421 m) summit

More information about the trail leading up to this landmark can be found at [https://utah.com/hiking/grandeur-peak](https://utah.com/hiking/grandeur-peak)

Grandeur Peak is a two part [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for sequencing of microbial isolates. Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

The "Grandeur" portion of this workflow is intended to be a species agnostic sequencing approach to paired-end Illumina sequencing quality control and assurance (QC) and serotyping in a local public health laboratory. "Grandeur" takes paired-end Illumina reads, removes adaptors and PHIX with [seqyclean](https://github.com/ibest/seqyclean), and creates contigs through _de novo_ alignment of the reads with [spades](https://cab.spbu.ru/software/spades/). Then a variety of QC and serotyping tools are employed. This workflow is similar to [Cutshaw](https://staph-b.github.io/staphb_toolkit/workflow_docs/cutshaw/), [Dryad](https://staph-b.github.io/staphb_toolkit/workflow_docs/dryad/), [Foushee](https://staph-b.github.io/staphb_toolkit/workflow_docs/foushee/), and [Tredegar](https://staph-b.github.io/staphb_toolkit/workflow_docs/tredegar/) in the [staphB toolkit](https://github.com/StaPH-B/staphb_toolkit) or [Bactopia](https://github.com/bactopia/bactopia), and the authors or maintainers of this workflow will not be offended if an alternative is used.

The "Peak" portion of this workflow uses the resulting contigs in fasta files (or in gff files created by [prokka](https://github.com/tseemann/prokka)) from the "Grandeur" portion through to a phylogenetic tree. The workflow, however, can be used with any fasta - including genomes downloaded from NCBI, or contigs generated from another workflow. In other words, the end user doesn't need to run "Grandeur" to run "Peak". The gff files created via [prokka](https://github.com/tseemann/prokka), are used by [roary](https://sanger-pathogens.github.io/Roary/) to define and align a core genome (a core genome is the genes that all the provided isolates share, also known as a [pan-genome](https://en.wikipedia.org/wiki/Pan-genome)). This multiple sequence alignment is used to create a tree with [iqtree2](http://www.iqtree.org/) and SNP matrix with [snp-dists](https://github.com/tseemann/snp-dists). As this workflow is dependent on the core genome, all of the fasta files put into this portion of the workflow must be the _same species_ (or predicted to have some sort of similar origin, so related plasmids will work as well). [@erinyoung](https://github.com/erinyoung) strongy recommends that the end user checks 'peak/roary/summary_statistics.txt' after running to ensure that the number of genes in the core genome makes sense. 

"Grandeur Peak" will also probably be a workflow of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit) once [@erinyoung](https://github.com/erinyoung) gets around to it and all the containers are ready.

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
- [Singularity](https://singularity.lbl.gov/install-linux) or [Docker](https://docs.docker.com/get-docker/) 

# Running Grandeur

![alt text](./configs/grandeur.png)

```
nextflow run grandeur.nf -c configs/singularity.config
```

"Grandeur" is expecting paired-end Illumina sequences of isolates at `workflow.launchDir + '/reads'` (`workflow.launchDir` is the directory where nextflow command is being run). If the directory of reads is at some other location, it can be specified with `params.reads` in a config file or on the command line as
```
nextflow run grandeur.nf -c configs/singularity.config --reads <directory of reads>
```
A directory will produce files at `workflow.launchDir + '/grandeur'`, but this can also be adjusted with `params.outdir` the same way.

"Grandeur" can also take fastas. Thus, a fasta can be created by a different workflow, like [Donut Falls](https://github.com/UPHL-BioNGS/Donut_Falls), and then go through this workflow for relevant serotyping and QC information. [Shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper) and [blobtools](https://blobtools.readme.io/docs) require fastq files, so these tools will not run on incoming fasta files. [Mash](https://github.com/marbl/Mash), [seqsero2](https://github.com/denglab/SeqSero2), [quast](http://quast.sourceforge.net/quast), [prokka](https://github.com/tseemann/prokka), [AMRfinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/), [serotypefinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/), [kleboarte](https://github.com/katholt/Kleborate), [kraken2](https://ccb.jhu.edu/software/kraken2/), and [mlst](https://github.com/tseemann/mlst) will work on fastas. The only parameters that may need to be configured by the **End User** is to specify the directory where the fasta files are with `'params.fasta'` (the default is `workflow.launchDir + '/fastas'`).
```
nextflow run grandeur.nf -c configs/singularity.config --fastas <directory of fastas>
```

This workflow will automatically grab any fastq files in `workflow.launchDir + '/reads'` and and fasta files in `workflow.launchDir + '/fastas'`. This is a feature of the workflow.

Additionally, "Grandeur" can use two optional tools to find contamination : [kraken2](https://ccb.jhu.edu/software/kraken2/) or [blobtools](https://blobtools.readme.io/docs) (or both!). These both use large databases that should be downloaded separately. 

## If the **End User** would like to use [kraken2](https://ccb.jhu.edu/software/kraken2/) to identify contamination:

First, there needs to be a kraken2 database. 

```
mkdir kraken2_db
cd kraken2_db
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz
tar -zxvf minikraken2_v2_8GB_201904.tgz
```
This downloads and expands the directory to 'kraken2_db/minikraken2_v2_8GB_201904_UPDATE'.

Then the corresponding params need to be updated.
- `params.kraken2` needs to be set to `true`
  - `params.kraken2 = 'true'`
- `params.kraken2_db` must be set to the directory with the kraken2 database
  - `params.kraken2_db = '<kraken2 database directory>'`
  - `params.kraken2_db = 'kraken2_db/minikraken2_v2_8GB_201904_UPDATE'` in the example

## UNDER DEVELOPMENT : If the End User would like to use [blobtools](https://blobtools.readme.io/docs) to identify contamination:

[Blobtools](https://blobtools.readme.io/docs) uses a blast database, so there needs to be a blast database:

```
mkdir blast_db
cd blast_db
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz"
for file in *tar.gz ; do tar -zxvf $file ; done
```
This downloads NCBI's NT database into 'blast_db'. 
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
├── prokka                                                                         # optional, but may save time by pre-generating gff files
│   └── sample
│       ├── sample.err
│       ├── sample.faa
│       ├── sample.ffn
│       ├── sample.fna
│       ├── sample.fsa
│       ├── sample.gbk
│       ├── sample.gff                                                             # annotated contig file that can be used via roary
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
├── shuffled
│   └── sample_shuffled.fastq.gz
├── spades
│   └── sample
│       ├── assembly_graph_after_simplification.gfa
│       ├── assembly_graph.fastg
│       ├── assembly_graph_with_scaffolds.gfa
│       ├── before_rr.fasta
│       ├── contigs.fasta
│       ├── contigs.paths
│       ├── dataset.info
│       ├── input_dataset.yaml
│       ├── K127
│       │   ├── assembly_graph_after_simplification.gfa
│       │   ├── assembly_graph.fastg
│       │   ├── assembly_graph_with_scaffolds.gfa
│       │   ├── before_rr.fasta
│       │   ├── configs
│       │   │   ├── careful_mda_mode.info
│       │   │   ├── careful_mode.info
│       │   │   ├── config.info
│       │   │   ├── construction.info
│       │   │   ├── detail_info_printer.info
│       │   │   ├── distance_estimation.info
│       │   │   ├── hmm_mode.info
│       │   │   ├── isolate_mode.info
│       │   │   ├── large_genome_mode.info
│       │   │   ├── mda_mode.info
│       │   │   ├── meta_mode.info
│       │   │   ├── metaplasmid_mode.info
│       │   │   ├── metaviral_mode.info
│       │   │   ├── moleculo_mode.info
│       │   │   ├── pe_params.info
│       │   │   ├── plasmid_mode.info
│       │   │   ├── rna_mode.info
│       │   │   ├── rnaviral_mode.info
│       │   │   ├── simplification.info
│       │   │   ├── toy.info
│       │   │   └── tsa.info
│       │   ├── final_contigs.fasta
│       │   ├── final_contigs.paths
│       │   ├── final.lib_data
│       │   ├── path_extend
│       │   ├── scaffolds.fasta
│       │   └── scaffolds.paths
│       ├── K21
│       │   ├── configs
│       │   │   ├── careful_mda_mode.info
│       │   │   ├── careful_mode.info
│       │   │   ├── config.info
│       │   │   ├── construction.info
│       │   │   ├── detail_info_printer.info
│       │   │   ├── distance_estimation.info
│       │   │   ├── hmm_mode.info
│       │   │   ├── isolate_mode.info
│       │   │   ├── large_genome_mode.info
│       │   │   ├── mda_mode.info
│       │   │   ├── meta_mode.info
│       │   │   ├── metaplasmid_mode.info
│       │   │   ├── metaviral_mode.info
│       │   │   ├── moleculo_mode.info
│       │   │   ├── pe_params.info
│       │   │   ├── plasmid_mode.info
│       │   │   ├── rna_mode.info
│       │   │   ├── rnaviral_mode.info
│       │   │   ├── simplification.info
│       │   │   ├── toy.info
│       │   │   └── tsa.info
│       │   ├── final.lib_data
│       │   └── simplified_contigs
│       │       ├── contigs_info
│       │       ├── contigs.off
│       │       └── contigs.seq
│       ├── K33
│       │   ├── configs
│       │   │   ├── careful_mda_mode.info
│       │   │   ├── careful_mode.info
│       │   │   ├── config.info
│       │   │   ├── construction.info
│       │   │   ├── detail_info_printer.info
│       │   │   ├── distance_estimation.info
│       │   │   ├── hmm_mode.info
│       │   │   ├── isolate_mode.info
│       │   │   ├── large_genome_mode.info
│       │   │   ├── mda_mode.info
│       │   │   ├── meta_mode.info
│       │   │   ├── metaplasmid_mode.info
│       │   │   ├── metaviral_mode.info
│       │   │   ├── moleculo_mode.info
│       │   │   ├── pe_params.info
│       │   │   ├── plasmid_mode.info
│       │   │   ├── rna_mode.info
│       │   │   ├── rnaviral_mode.info
│       │   │   ├── simplification.info
│       │   │   ├── toy.info
│       │   │   └── tsa.info
│       │   ├── final.lib_data
│       │   └── simplified_contigs
│       │       ├── contigs_info
│       │       ├── contigs.off
│       │       └── contigs.seq
│       ├── K55
│       │   ├── configs
│       │   │   ├── careful_mda_mode.info
│       │   │   ├── careful_mode.info
│       │   │   ├── config.info
│       │   │   ├── construction.info
│       │   │   ├── detail_info_printer.info
│       │   │   ├── distance_estimation.info
│       │   │   ├── hmm_mode.info
│       │   │   ├── isolate_mode.info
│       │   │   ├── large_genome_mode.info
│       │   │   ├── mda_mode.info
│       │   │   ├── meta_mode.info
│       │   │   ├── metaplasmid_mode.info
│       │   │   ├── metaviral_mode.info
│       │   │   ├── moleculo_mode.info
│       │   │   ├── pe_params.info
│       │   │   ├── plasmid_mode.info
│       │   │   ├── rna_mode.info
│       │   │   ├── rnaviral_mode.info
│       │   │   ├── simplification.info
│       │   │   ├── toy.info
│       │   │   └── tsa.info
│       │   ├── final.lib_data
│       │   └── simplified_contigs
│       │       ├── contigs_info
│       │       ├── contigs.off
│       │       └── contigs.seq
│       ├── K77
│       │   ├── configs
│       │   │   ├── careful_mda_mode.info
│       │   │   ├── careful_mode.info
│       │   │   ├── config.info
│       │   │   ├── construction.info
│       │   │   ├── detail_info_printer.info
│       │   │   ├── distance_estimation.info
│       │   │   ├── hmm_mode.info
│       │   │   ├── isolate_mode.info
│       │   │   ├── large_genome_mode.info
│       │   │   ├── mda_mode.info
│       │   │   ├── meta_mode.info
│       │   │   ├── metaplasmid_mode.info
│       │   │   ├── metaviral_mode.info
│       │   │   ├── moleculo_mode.info
│       │   │   ├── pe_params.info
│       │   │   ├── plasmid_mode.info
│       │   │   ├── rna_mode.info
│       │   │   ├── rnaviral_mode.info
│       │   │   ├── simplification.info
│       │   │   ├── toy.info
│       │   │   └── tsa.info
│       │   ├── final.lib_data
│       │   └── simplified_contigs
│       │       ├── contigs_info
│       │       ├── contigs.off
│       │       └── contigs.seq
│       ├── K99
│       │   ├── configs
│       │   │   ├── careful_mda_mode.info
│       │   │   ├── careful_mode.info
│       │   │   ├── config.info
│       │   │   ├── construction.info
│       │   │   ├── detail_info_printer.info
│       │   │   ├── distance_estimation.info
│       │   │   ├── hmm_mode.info
│       │   │   ├── isolate_mode.info
│       │   │   ├── large_genome_mode.info
│       │   │   ├── mda_mode.info
│       │   │   ├── meta_mode.info
│       │   │   ├── metaplasmid_mode.info
│       │   │   ├── metaviral_mode.info
│       │   │   ├── moleculo_mode.info
│       │   │   ├── pe_params.info
│       │   │   ├── plasmid_mode.info
│       │   │   ├── rna_mode.info
│       │   │   ├── rnaviral_mode.info
│       │   │   ├── simplification.info
│       │   │   ├── toy.info
│       │   │   └── tsa.info
│       │   ├── final.lib_data
│       │   └── simplified_contigs
│       │       ├── contigs_info
│       │       ├── contigs.off
│       │       └── contigs.seq
│       ├── misc
│       │   └── broken_scaffolds.fasta
│       ├── params.txt
│       ├── pipeline_state
│       │   ├── stage_0_before_start
│       │   ├── stage_10_bs
│       │   ├── stage_11_terminate
│       │   ├── stage_1_as_start
│       │   ├── stage_2_k21
│       │   ├── stage_3_k33
│       │   ├── stage_4_k55
│       │   ├── stage_5_k77
│       │   ├── stage_6_k99
│       │   ├── stage_7_k127
│       │   ├── stage_8_copy_files
│       │   └── stage_9_as_finish
│       ├── run_spades.sh
│       ├── run_spades.yaml
│       ├── scaffolds.fasta
│       ├── scaffolds.paths
│       ├── spades.log
│       └── tmp
└── summary
    ├── sample.summary.tsv
    ├── sample.summary.txt
    └── grandeur_summary.txt                                                    # a table with a summary from all the serotyping and QC tools
```

</details>

# Running Peak

![alt text](./configs/Peak.png)

```
nextflow run peak.nf -c configs/singularity.config
```

"Peak" is expecting fasta files, [prokka](https://github.com/tseemann/prokka)'s gff files (which are fasta files that have been annotated with [prokka](https://github.com/tseemann/prokka)), or both. Any file that hasn't already gone through [prokka](https://github.com/tseemann/prokka) will be annotated by [prokka](https://github.com/tseemann/prokka). Then the gff files from [prokka](https://github.com/tseemann/prokka) are aligned with [roary](https://sanger-pathogens.github.io/Roary/). The multiple sequence alignment of the core genome is put through [iqtree2](http://www.iqtree.org/) and [snp-dists](https://github.com/tseemann/snp-dists). 

Fasta files are expected at `workflow.launchDir + '/fastas'`, and can be adjusted with `'params.fastas'`.
GFF files are expected `workflow.launchDir + '/gff'`, and can be adjusted with `'params.gff'`. 

On the command line it looks something like this:
```
nextflow run peak.nf -c configs/singularity.config --fastas <directory of fasta files> --gff <directory of gff files>
```

## Although not required, [@erinyoung](https://github.com/erinyoung) suggests changing some of the default parameters. 
- Changing the prokka parameters to include the genus and species of your organism as well as your 'centre'.  
  - `params.prokka_options` = `'--genus <genus> --species <species>'`
  - `params.center` = `'<center>'` 
- Changing the iqtree2 parameters to designate which file should be used as an outgroup
  - `params.outgroup` - `'<outgroup>'`

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
- [spades](https://cab.spbu.ru/software/spades/) - _de novo_ alignment
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
There is a template file with all the variables in this repo at [configs/grandeur_template.config](./configs/grandeur_template.config) that the **End User** can copy and edit. All of the parameters are included in that file.

There's also a config file what we use here at UPHL, [UPHL.config](./configs/UPHL.config).

## Do you have test data?
Yes, actually, but also not quite. Fastq files can be rather large and github isn't the best place to store them. Instead, some sample fastq files from the [SRA](https://www.ncbi.nlm.nih.gov/sra) recommended. These have gone through the "Grandeur" workflow and the results are included in this repo for your comparison. 
- [SRR7889058](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7889058) : _Salmonella enterica_
- [SRR11725329](https://www.ncbi.nlm.nih.gov/sra/?term=SRR11725329) : _Escherichia coli_
- [SRR7738178](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7738178) : _Shigella sonnei_ 
- [SRR14634837](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14634837) : _Klebsiella pneumoniae_
- [SRR13643280](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13643280) : _Acinetobacter baumannii_

These can be downloaded using [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) (part of the [SRA](https://github.com/ncbi/sra-tools) toolkit) as follows:
```
fasterq-dump -A SRR7889058 --split-files --outdir reads
```

Followed by
```
nextflow grandeur.nf -c configs/singularity.config
```

The directory `data/fasta` contains fasta files for 
- _Stenotrophomonas maltophilia_ (GCF_900475405.1_44087_C01_genomic)
- _Acinetobacter baumannii_ (GCF_008632635.1_ASM863263v1_genomic)
- _Klebsiella pneumoniae_ (GCF_000240185.1_ASM24018v2_genomic)
- _Shigella flexneri_ (GCF_000006925.2_ASM692v2_genomic)
- _Salmonella enterica_ (GCF_000006945.2_ASM694v2_genomic)
- _Escherichia coli_ (GCF_000005845.2_ASM584v2_genomic)

The workflow can be tested on these fastas files using the following command.

```
nextflow grandeur.nf -c configs/singularity.config --fastas data/fasta
```
Summary files from running these through the default workflow ([grandeur_results.tsv](./data/grandeur_results.tsv)) as well as with UPHL's config file ([UPHL_grandeur_results.tsv](./data/UPHL_grandeur_results.tsv)) are also available.


The directory `data/peak_test` contains one gff file and 6 fasta files of _Stenotrophomonas maltophilia_ that can be used to test "Peak". A resulting treefile ([iqtree.treefile](./data/peak_test/iqtree.treefile)), snp_matrix ([snp_matrix.txt](./data/peak_test/snp_matrix.txt)), and roary summary file ([summary_statistics.txt](./data/peak_test/summary_statistics.txt)) are included for comparison. 

The basics of testing out "Peak":
```
nextflow peak.nf -c configs/singularity.config --fastas data/peak_test --gff data/peak_test
```

## What about CLIA validation?
At UPHL, we use this workflow to determine the serotype of Salmonella and E. coli under CLIA. Therefore, if you look at [our config file](./configs/UPHL.config), we explicitly specify the containers used for Salmonella serotyping with [seqsero2](https://github.com/denglab/SeqSero2), as well as E. coli serotyping with [serotypefinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/) and [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper). We also specify the containers used prior to these processes in the workflow : [seqyclean](https://github.com/ibest/seqyclean) and [spades](https://cab.spbu.ru/software/spades/) (which isn't used for [seqsero2](https://github.com/denglab/SeqSero2) or [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper), but it is used for [serotypefinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/)). 

The CLIA officer of the **End User** may request additional locks be put in place, like having all of the containers specified. If additional help is needed, please [submit an issue](https://github.com/UPHL-BioNGS/Grandeur/issues) or [Email me](eriny@utah.gov).

## How were serotyping tools chosen for this workflow?
They perform _well_, their containers were easy to create, and [@erinyoung](https://github.com/erinyoung) had heard about them.

## Are any other tools getting added to "Grandeur"?

As "Grandeur" is intended to be a species agnostic workflow for a local public health laboratory, and sequencing is continuing to expand in its utility, new tools are constantly being needed to analyze isolates to further public health goals. 

For example, the _C. auris_ worfklow [mycosnp](https://github.com/CDCgov/mycosnp) is missing from this workflow, as is any sort of plasmid analysis. 

Many of these additional tools are added by need locally or from the **End User**, so if the **End User** knows of other serotyping/analysis tools, please [submit an issue](https://github.com/UPHL-BioNGS/Grandeur/issues) or tell [@erinyoung](https://github.com/erinyoung) about it, and we'll work in some options. 

[@erinyoung](https://github.com/erinyoung) also appreciates pull requests from forks.

**Warning** : If there's not a relaible container of the suggested tool, [@erinyoung](https://github.com/erinyoung) will request that the **End User** create a container for that tool and contribute to [StaPH-B's docker repositories](https://github.com/StaPH-B/docker-builds).

## What about organisms with large genomes?
Organisms with large genomes can still contribute to disease, but this is not the workflow for those. "Grandeur" uses [spades](https://cab.spbu.ru/software/spades/). Large genomes may be too much for [spades](https://cab.spbu.ru/software/spades/). 

## What about SARS-CoV-2?
As of the time of writing this README, reference-based alignment of SARS-CoV-2 is still the norm. "Grandeur" is for _de novo_ assembly of things with small genome. [Cecret](https://github.com/UPHL-BioNGS/Cecret) would be a better workflow for SARS-CoV-2 sequencing. 

## What is genome_sizes.json used for?
[genome_sizes.json](./configs/genome_sizes.json) has a list of commonly sequenced organisms and the approximate expected genome size for each organism. This is only used for the "cg-pipeline" process to estimate coverage. A file from the **End User** can be used instead and specified with `params.genome_sizes`. 

## How do I cite this workflow?

This workflow stands on the shoulders of giants. As such, please cite the individual tools that were useful for your manuscript so that those developers can continue to get funding. They are listed above. Mentioning this workflow in the text as "The Grandeur workflow v.VERSION (www.github.com/UPHL-BioNGS/Grandeur)" is good enough for [@erinyoung](https://github.com/erinyoung)'s ego.

## Why doesn't the default [shigatyper](https://hub.docker.com/r/andrewlangvt/shigatyper) container work with nextflow's `-with-tower` option?

There is actually hope that this will eventually be remedied, but for now the default shigatyper container does not contain 'ps', so it will not work with nextflow tower. 

## Can I use roary's QC options with kraken?

Yes. There are params that can be adjusted to use a kraken database for roary's qc options. Adjust `params.kraken = true` and `params.kraken_db = 'kraken_db'`. This does not work with staphb's current [roary](https://hub.docker.com/r/staphb/roary) container.

![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)


