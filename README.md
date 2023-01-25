# Grandeur

<img src="https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg" width="500" align="left" />

Named after the beautiful [Grandeur Peak](https://www.alltrails.com/trail/us/utah/grandeur-peak-east-trail-from-church-fork)

Image Credit: [ryancornia](https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg)

Location:  40.707, -111.76, 8,299 ft (2,421 m) summit

More information about the trail leading up to this landmark can be found at [https://utah.com/hiking/grandeur-peak](https://utah.com/hiking/grandeur-peak)

Grandeur is a [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/). "Grandeur" is intended to be a species agnostic sequencing analysis workflow to paired-end Illumina sequencing quality control and assurance (QC) and serotyping in a local public health laboratory. 

"Grandeur" is meant to augment CDC's [PHOENIX](https://github.com/CDCgov/phoenix) nextflow workflow, which is the recommended usage. In principle, the cleaned reads and contigs from PHOENIX undergo additional quality metric and serotyping steps, with a heavy emphasis on [fastANI](https://github.com/ParBLiSS/FastANI) and [AMRFinderPlus](https://github.com/ncbi/amr).

"Grandeur" can be a standalone workflow that takes paired-end Illumina reads, removes adaptors with [fastp](https://github.com/OpenGene/fastp) and PHIX with [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), and creates contigs through _de novo_ alignment of the reads with [spades](https://cab.spbu.ru/software/spades/).

"Grandeur" can also take in fastas. Additional fastas includes genomes downloaded from NCBI or contigs generated from another workflow, such as [Donut Falls](https://github.com/UPHL-BioNGS/Donut_Falls). The gff files created via [prokka](https://github.com/tseemann/prokka), are used by [roary](https://sanger-pathogens.github.io/Roary/) to define and align a core genome (a core genome is the genes that all the provided isolates share, also known as a [pan-genome](https://en.wikipedia.org/wiki/Pan-genome)). This multiple sequence alignment is used to create a tree with [iqtree2](http://www.iqtree.org/) and SNP matrix with [snp-dists](https://github.com/tseemann/snp-dists). As this workflow is dependent on the core genome, all of the fasta files put into this sub-workflow must be the _same species_ (or predicted to have some sort of similar origin, so related plasmids will work as well). [@erinyoung](https://github.com/erinyoung) strongy recommends that the end user checks `'grandeur/roary/summary_statistics.txt'` after running to ensure that the number of genes in the core genome makes sense.

"Grandeur" is also a workflow of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit)

## Dependencies

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Singularity](https://singularity.lbl.gov/install-linux) or [Docker](https://docs.docker.com/get-docker/)
- [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

## Usage

Default workflow that takes fastq files, runs them through QC/serotyping/etc, creates contig files
```
# using singularity
nextflow run UPHL-BioNGS/Grandeur -profile singularity --reads <path to reads>
# using docker
nextflow run UPHL-BioNGS/Grandeur -profile docker --reads <path to reads>
```

### Choosing the right profile
For simplificity, Grandeur has some profiles which _should_ fit the majority of uses. The workflow is meant to work with containers and has basic profiles for both docker and singularity, which are in the example above. Otherwise, *ALL* the commands used in the workflow will need to be present PATH. 

#### Choose a container manager
- singularity : use singularity to manage containers
- docker : use docker to manage containers

### Choose the start and end point
- msa                    : for multiple sequence alignment with roary (all inputs should be related) of input files
- just_msa               : for multiple sequence alignment with roary (all inputs should be related) of input files, and turns off processes not directly used
- uphl                   : the profile used at UPHL (is not intended to work on other systems)

WARNING: All input files for `*msa*` profiles must all be somewhat related (i.e. same species) because they need to share enough genes in their core genome.

## Default file structure
(can be adjusted with 'params.reads' and 'params.fastas')

### Paired-end fastq.gz (ending with 'fastq', 'fastq.gz', 'fq', or 'fq.gz') reads as follows or designate directory with 'params.reads' or '--reads'
```
directory
└── reads
     └── *fastq.gz
```
WARNING : Sometimes nextflow does not catch every name of paired-end fastq files. This workflow is meant to be fairly agnostic, but if paired-end fastq files are not being found it might be worth renaming them to some sort of 'sample_1.fastq.gz', 'sample_2.fastq.gz' format.

### Fasta files (ending with 'fa', 'fasta', or 'fna') as follows or designate directory with 'params.fastas' or '--fastas'
```
directory
└── fastas
     └── *fasta
```

WARNING : "Grandeur" will automatically grab any fastq files in `workflow.launchDir + '/reads'`, fasta files in `workflow.launchDir + '/fastas'`, and gff files in `workflow.launchDir + '/gff'`. This is a **feature** of the workflow.

## The directory where the results will be located
A directory will produce files at `'grandeur'` in where the command was inputted, but this can also be adjusted with 'params.outdir' or '--outdir'.

<details>
   <summary>Final file tree after running grandeur.nf</summary>

```
grandeur/
├── aligned
│   ├── sample.sorted.bam
│   └── sample.sorted.bam.csi
├── bbduk
│   ├── sample.matched_phix.fq
│   ├── sample.phix.stats.txt
│   ├── sample_rmphix_R1.fastq.gz
│   └── sample_rmphix_R2.fastq.gz
├── blastn
│   └── sample.tsv
├── blobtools
│   ├── sample.sample.sorted.bam.cov
│   ├── sample.blobDB.json
│   ├── sample.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png
│   ├── sample.blobDB.json.bestsum.species.p8.span.100.blobplot.read_cov.bam0.png
│   ├── sample.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt          # Genus and species of the reads
│   └── sample.blobDB.table.txt
├── contigs
│   └── sample_contigs.fa                                                          # fasta file of contigs
├── sample
│   ├── fastani.out
│   └── sample.txt
├── fastp
│   ├── sample_fastp.html
│   ├── sample_fastp.json
│   ├── sample_fastp_R1.fastq.gz
│   └── sample_fastp_R2.fastq.gz
├── fastqc
│   ├── sample_fastqc.html
│   ├── sample_fastqc.zip
│   ├── sample_fastqc.html
│   └── sample_fastqc.zip
├── gff
│   └── sample.gff                                                                 # gff file created by prokka
├── grandeur_results.tsv                                                           # summary file
├── iqtree2
│   ├── iqtree.ckp.gz
│   ├── iqtree.contree                                                             # treefile without node values
│   ├── iqtree.iqtree
│   ├── iqtree.log
│   ├── iqtree.splits.nex
│   └── iqtree.treefile
├── kleborate
│   ├── sample_results.txt
│   └── kleborate_results.txt                                                      # klebsiella hypervirulence scoring
├── kraken2
│   └── sample_kraken2_report.txt
├── logs
├── mash
│   ├── sample_mashdist.txt                                                        # mash distances
│   └── sample.msh
├── mlst
│   ├── sample_mlst.txt
│   └── mlst_result.tsv                                                            # mlst of organism (if found)
├── multiqc
│   ├── multiqc_data
│   │   └── *
│   └── multiqc_report.html
├── ncbi-AMRFinderplus
├── plasmidfinder
│   └── sample
│       ├── data.json
│       └── tmp
│           ├── out_enterobacteriaceae.xml
│           ├── out_Inc18.xml
│           ├── out_NT_Rep.xml
│           ├── out_Rep1.xml
│           ├── out_Rep2.xml
│           ├── out_Rep3.xml
│           ├── out_RepA_N.xml
│           ├── out_RepL.xml
│           └── out_Rep_trans.xml
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
│   ├── core_gene_alignment.aln                                                 # core genome alignment
│   ├── fixed_input_files
│   │   └── sample.gff
│   ├── gene_presence_absence.csv
│   ├── gene_presence_absence.Rtab
│   ├── number_of_conserved_genes.Rtab
│   ├── number_of_genes_in_pan_genome.Rtab
│   ├── number_of_new_genes.Rtab
│   ├── number_of_unique_genes.Rtab
│   ├── pan_genome_reference.fa
│   └── summary_statistics.txt                                                 # important file with the number of genes involved in core genome
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
├── snp_dists
│   └── snp_matrix.txt                      # SNP matrix counting the number of SNPs that each sample differs by
└── summary
    ├── sample.summary.tsv
    ├── sample.summary.txt
    └── grandeur_summary.txt                # a table with a summary from all the serotyping and QC tools
```
</details>

## Finding contamination
Additionally, "Grandeur" can use two optional tools to find cross-species contamination : [kraken2](https://ccb.jhu.edu/software/kraken2/) or [blobtools](https://blobtools.readme.io/docs) (or both!). These both use large databases that should be downloaded separately.

### If the **End User** would like to use [kraken2](https://ccb.jhu.edu/software/kraken2/) to identify contamination (optional):

First, there needs to be a kraken2 database.

```
mkdir kraken2_db
cd kraken2_db
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz
tar -zxvf minikraken2_v2_8GB_201904.tgz
```
This downloads and expands the directory to 'kraken2_db/minikraken2_v2_8GB_201904_UPDATE'.

Then the corresponding kraken2_db param needs to be updated. `params.kraken2_db` must be set to the directory with the kraken2 database. The config file lines for the above example (can be copied and pasted into a config file):
```
params.kraken2_db = 'kraken2_db/minikraken2_v2_8GB_201904_UPDATE'
```

### If the **End User** would like to use [blobtools](https://blobtools.readme.io/docs) to identify contamination (optional):

[Blobtools](https://blobtools.readme.io/docs) uses a blast database, so there needs to be a blast database:

```
mkdir blast_db
cd blast_db

# get the taxdump file
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
# get the nt files
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz"
# decompress the nt files
for file in *tar.gz ; do tar -zxvf $file ; done
```
This downloads NCBI's NT database into 'blast_db'.
WARNING: The download and decompressing can take a long time.

- `params.blast_db` must be set to the directory with the blast database
- `params.blast_db_type` must be set to the type of blast database that is being used

The config file lines for the above example (can be copied and pasted into a config file):
```
params.blast_db = 'blast_db'
params.blast_db_type = 'nt'
```

* at UPHL, we use the 'ref_prok_rep_genomes' blast database instead of 'nt'

## Updating the mash reference file (optional)

The mash reference file in the staphb/mash container is from RefSeq77. The most current RefSeq version for January 11, 2023 (the time of updating this readme) is 215. Many organism names and even some species defintions have changed. As such, the **End User** can create a custom, more-current mash database to use with this workflow. UPHL has created an '*.msh' from RefSeq version 215 and is available for download at []().

```
wget URL
```

The config file lines for the above example (can be copied and pasted into a config file):
```
params.mash_db = "rep-genomes.msh"
```


![alt text](./configs/grandeur_lucidchart.png)


## Suggested inputs

Although not required, [@erinyoung](https://github.com/erinyoung) suggests changing `params.iqtree2_outgroup` when creating a phylogenetic tree with iqtree2. For a file with a name like 'GCF_000006765.1_ASM676v1_genomic.gff', the outgroup would be 'GCF_000006765.1_ASM676v1_genomic'.

Example lines for a config file:
```
iqtree2_outgroup = 'GCF_000006765.1_ASM676v1_genomic'
```

## Examples of command line usage
Phew! With all of that out of the way, what does it look like on the command line?

```
# vanilla (uses refseq 77 with mash, downloads genomes from NCBI)
nextflow run UPHL-BioNGS/Grandeur -profile singularity --read <path to reads>

# without downloading genomes from NCBI (uses refseq 77 with mash)
nextflow run UPHL-BioNGS/Grandeur -profile singularity --read reads --current_datasets false

# running on a collection of fastq files and fasta files that will be placed in a phylogenetic tree
# (will download genomes from NCBI and add the FastANI top hit to tree building)
nextflow run UPHL-BioNGS/Grandeur -profile singularity,msa --reads reads --fastas fastas

# only wanting to create a tree for some fastas
nextflow run UPHL-BioNGS/Grandeur -profile singularity,just_msa --fasta fastas --iqtree2_outgroup outgroup
```

# Visualizing the tree


It would be nice if an image of tree was automatically generated from this. Really nice. It has been difficult to find a container that creates a high-quality tree in a command line interface that automatically resizes the text and tree appropriately. (Although, we're always open to suggestions!) The phylogenetic tree found at 'grandeur/iqtree2/iqtree.treefile' or 'grandeur/iqtree2/iqtree.contree' can be visualized through multiple tools, such as [ggtree](https://yulab-smu.top/treedata-book/), [microreact](https://microreact.org/), or [itol](https://itol.embl.de/).  


# Grandeur wouldn't be possible without:

- [amrfinderplus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) - identification of genes associated with antimicrobial resistence
- [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) - removal of PhiX
- [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) - maps reads to contigs for blobtools
- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK279684/) - read identification with blobtools
- [blobtools](https://blobtools.readme.io/docs) - contamination
- [datasets](https://github.com/ncbi/datasets) - downloads genomes from NCBI
- [fastani](https://github.com/ParBLiSS/FastANI) - species evaluator
- [fastp](https://github.com/OpenGene/fastp) - cleaning reads
- [fastqc](https://github.com/s-andrews/FastQC) - fastq file QC
- [fastq-scan](https://github.com/rpetit3/fastq-scan) - quality metrics for fastq files in json format
- [iqtree2](http://www.iqtree.org/) - phylogenetic tree creation - used after core genome alignment
- [kleborate](https://github.com/katholt/Kleborate) - Klebsiella serotyping
- [kraken2](https://ccb.jhu.edu/software/kraken2/) - contamination
- [mash](https://github.com/marbl/Mash) - species identifier
- [mlst](https://github.com/tseemann/mlst) - identification of MLST subtype
- [multiqc](https://multiqc.info/) - summarizes QC efforts
- [plasmidfinder](https://bitbucket.org/genomicepidemiology/plasmidfinder/) - MLST typing for plasmids
- [prokka](https://github.com/tseemann/prokka) - gene annotation - used for core genome alignment
- [quast](http://quast.sourceforge.net/quast) - contig QC
- [roary](https://sanger-pathogens.github.io/Roary/) - core genome alignment - optional
- [seqsero2](https://github.com/denglab/SeqSero2) - Salmonella serotyping
- [serotypefinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/) - E. coli serotyping
- [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper) - Shigella serotyping
- [snp-dists](https://github.com/tseemann/snp-dists) - SNP matrix - used after core genome aligment
- [spades](https://cab.spbu.ru/software/spades/) - _de novo_ alignment
- [sra toolkit](https://hpc.nih.gov/apps/sratoolkit.html) - for downloading fastq files for testing


