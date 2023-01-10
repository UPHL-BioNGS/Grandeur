#!/usr/bin/env nextflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Welcome to this workflow! Issues and contributions are gladly accepted at https://github.com/UPHL-BioNGS/Grandeur .

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

println("Currently using the Grandeur workflow for use with microbial sequencing. The view is great from 8299 feet (2530 meters) above sea level.\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: ${workflow.manifest.version}")
println("")

nextflow.enable.dsl               = 2

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Defining params

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

params.outdir                     = 'grandeur'
params.maxcpus                    = 12
params.medcpus                    = 4
params.minimum_reads              = 10000

// connecting to phoenix
params.phoenix_wf                 = false

// connecting to Donut Falls
params.donut_falls_wf             = false

// input files
params.reads                      = workflow.launchDir + '/reads'
params.fastas                     = workflow.launchDir + '/fastas'
params.gff                        = workflow.launchDir + '/gff'

// external files
params.kraken2_db                 = ''
params.blast_db                   = ''
params.mash_db                    = ''
params.plasmidfinder_ref          = ''
params.fastani_ref                = workflow.projectDir + "/configs/fastani_ref.tar.gz"
params.genome_sizes               = workflow.projectDir + "/configs/genome_sizes.json"

// for testing
params.sra_accessions             = []

// tool-specific command line options
params.amrfinderplus_options      = ''
params.bbduk_options              = 'k=31 hdist=1'
params.bbmap_options              = ''
params.blast_db_type              = 'nt'
params.blastn_options             = '-max_target_seqs 10 -max_hsps 1 -evalue 1e-25'
params.blobtools_create_options   = ''
params.blobtools_view_options     = ''
params.blobtools_plot_options     = '--format png -r species'
params.blobtools_bbmap_options    = ''
params.cg_pipeline_options        = '--qual_offset 33 --minLength 1'
params.current_datasets           = true
params.datasets_max_genomes       = 5
params.fastani_options            = '--matrix'
params.fasterqdump_options        = ''
params.fastp_options              = '--detect_adapter_for_pe'
params.fastqc_options             = ''
params.fastqscan_options          = ''
params.iqtree2_options            = '-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000'
params.iqtree2_outgroup           = ''
params.kleborate_options          = '-all'
params.kraken2_options            = ''
params.mash_sketch_options        = '-m 2'
params.mash_dist_options          = '-v 0 -d 0.5'
params.mash_max_hits              = 100
params.mlst_options               = ''
params.multiqc_options            = ''
params.plasmidfinder_options      = ''
params.prokka_options             = '--mincontiglen 500 --compliant --locustag locus_tag --centre STAPHB'
params.quast_options              = ''
params.roary_options              = ''
params.roary_min_genes            = 1500
params.seqsero2_options           = '-m a -b mem'
params.serotypefinder_options     = ''
params.shigatyper_options         = ''
params.snp_dists_options          = '-c'
params.spades_options             = '--isolate'

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Sharing params with subworkflows

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

include { average_nucleotide_identity } from './subworkflows/average_nucleotide_identity.nf'  addParams(params)
include { blobtools }                   from './subworkflows/blobtools.nf'                    addParams(params)
include { de_novo_alignment }           from './subworkflows/de_novo_alignment.nf'            addParams(params)
include { information }                 from './subworkflows/information.nf'                  addParams(params)
include { kraken2 }                     from './subworkflows/kraken2.nf'                      addParams(params)
include { min_hash_distance }           from './subworkflows/min_hash_distance.nf'            addParams(params)
include { phylogenetic_analysis }       from './subworkflows/phylogenetic_analysis.nf'        addParams(params)
include { report }                      from './subworkflows/report.nf'                       addParams(params)
include { test }                        from './subworkflows/test.nf'                         addParams(params)

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for scripts

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Creating the summary files
summary_script = Channel.fromPath(workflow.projectDir + "/bin/summary.py", type:'file')

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for input files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Getting the fastq files
Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                  "${params.reads}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .view { "Paired-end fastq files found : ${it[0]}" }
  .unique()
  .set { ch_reads }

// Getting contig or fasta files
Channel
  .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
  .map { file -> tuple(file.baseName, file) }
  .unique()
  .set { ch_fastas }

// Getting fasta files that have been annotated with prokka
Channel.fromPath("${params.gff}/*.gff", type: 'file')
  .view { "gff file : $it" }
  .unique()
  .set { ch_gffs }

// Getting accession for downloading
ch_sra_accessions = Channel.from( params.sra_accessions )

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for database files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Getting the file with genome sizes of common organisms for cg-pipeline. The End User can use their own file and set with a param
ch_genome_sizes    = Channel.fromPath(params.genome_sizes, type:'file')

// Getting the reference genomes for fastANI
ch_fastani_genomes = Channel.fromPath("${params.fastani_ref}", type:'file')

// Getting the database for blobtools
ch_blast_db = params.blast_db
              ? Channel
                  .fromPath(params.blast_db, type:'dir')
                  .ifEmpty{
                    println("No blast database was found at ${params.blast_db}")
                    println("Set 'params.blast_db' to directory with blast database")
                    exit 1
                  }
                  .view { "Local Blast Database for Blobtools : $it" }
              : Channel.empty()

// Getting the kraken2 database
ch_kraken2_db = params.kraken2_db
              ? Channel
                  .fromPath(params.kraken2_db, type:'dir')
                  .ifEmpty{
                    println("No kraken2 database was found at ${params.kraken2_db}")
                    println("Set 'params.kraken2_db' to directory with kraken2 database")
                    exit 1
                  }
                  .view { "Local kraken2 database : $it" }
              : Channel.empty()

ch_mash_db = params.mash_db 
            ? Channel
              .fromPath(params.mash_db, type: 'file')
              .ifEmpty{
                println("No mash database was found at ${params.mash_db}")
                println("Set 'params.mash_db' to file of pre-sketched mash reference")
                exit 1
              }
              .view { "Mash reference : $it" }
            : Channel.empty()

println("The files and directory for results is " + params.outdir)
println("The maximum number of CPUS for any one process is ${params.maxcpus}")

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Workflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


workflow {

  ch_for_multiqc   = Channel.empty()
  ch_for_summary   = Channel.empty()
  ch_for_species   = Channel.empty()

  // getting test files
  if ( ! params.sra_accessions.isEmpty() ) { 
    test(ch_sra_accessions)
    ch_raw_reads = ch_reads.mix(test.out.fastq)
  } else {
    ch_raw_reads = ch_reads
  }

  if ( params.phoenix_wf )      { phoenix(input) }
  if ( params.donut_falls_wf )  { donut_falls(input) }

  // either phoenix or de_novo_alignment is required
  de_novo_alignment(ch_reads)
  
  ch_contigs       = ch_fastas.mix(de_novo_alignment.out.contigs)
  ch_clean_reads   = de_novo_alignment.out.clean_reads

  // optional subworkflow blobtools (useful for interspecies contamination)
  if (params.blast_db) {
    blobtools(ch_clean_reads, ch_contigs, ch_blast_db)
    ch_for_multiqc = ch_for_multiqc.mix(blobtools.out.for_multiqc)
    ch_for_summary = ch_for_summary.mix(blobtools.out.for_summary)
    ch_for_species = ch_for_species.mix(blobtools.out.for_species)
  }

  // optional subworkflow kraken2 (useful for interspecies contamination)
  if (params.kraken2_db) {
    kraken2(ch_clean_reads, ch_contigs, ch_kraken2_db)
    ch_for_multiqc = ch_for_multiqc.mix(kraken2.out.for_multiqc)
    ch_for_summary = ch_for_summary.mix(kraken2.out.for_summary)
    ch_for_species = ch_for_species.mix(kraken2.out.for_species)
  } 

  // subworkflow mash for species determination
  min_hash_distance(ch_clean_reads, ch_contigs, ch_mash_db)

  // determining organisms in sample
  average_nucleotide_identity(
    ch_for_summary.collect(),
    ch_contigs,
    ch_fastani_genomes
  )
  ch_for_species   = ch_for_species.mix(min_hash_distance.out.for_species).mix(average_nucleotide_identity.out.for_species)
  
  // getting all the other information
  ch_for_size      = average_nucleotide_identity.out.for_size.join(min_hash_distance.out.for_size, by: 0)
  information(ch_reads, 
    ch_contigs, 
    ch_for_species, 
    ch_for_size.combine(ch_genome_sizes))
  ch_for_multiqc   = ch_for_multiqc.mix(information.out.for_multiqc)

  // optional subworkflow for comparing shared genes
  // phylogenetic_analysis(contigs, mash_species, mash_genus, gffs )

  // getting a summary of everything
  ch_for_multiqc   = ch_for_multiqc.mix(de_novo_alignment.out.for_multiqc)

  ch_for_summary
    .mix(min_hash_distance.out.for_summary)
    .mix(information.out.for_summary)
    .mix(average_nucleotide_identity.out.for_summary)
    .concat(summary_script)
    .flatten()
    .set { for_summary }

  report(
    ch_reads, 
    ch_fastas, 
    ch_for_multiqc.collect(), 
    for_summary.collect(), 
    de_novo_alignment.out.fastp_reads, 
    de_novo_alignment.out.phix_reads)
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Final Steps

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("MultiQC report can be found at ${params.outdir}/multiqc/multiqc_report.html")
    println("Summary can be found at ${params.outdir}/grandeur_summary.tsv")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
