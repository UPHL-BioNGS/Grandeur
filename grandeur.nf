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

// Getting config file

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

params.config_file                = false
if (params.config_file) {
  def src = new File("${workflow.projectDir}/configs/grandeur_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println("A config file can be found at ${workflow.launchDir}/edit_me.config")
  exit 0
}

params.fastqscan_options = ""
if (params.fastqscan_options) {
  println("WARNING : ${params.fastqscan_options} is depreciated" )
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Defining params

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

params.outdir                     = "grandeur"
params.maxcpus                    = 12
params.medcpus                    = 4
params.minimum_reads              = 10000

// connecting to phoenix (in development)
// params.phoenix_wf                 = false
// params.phoenix_dir                = workflow.projectDir + "/../phoenix/workflows/phoenix.nf"
// params.input                      = ""

// connecting to Donut Falls (in development)
params.donut_falls_wf             = false

// input files
params.reads                      = workflow.launchDir + "/reads"
params.fastas                     = workflow.launchDir + "/fastas"
params.gff                        = workflow.launchDir + "/gff"
params.sample_sheet               = ""

// external files
params.kraken2_db                 = ""
params.blast_db                   = ""
params.mash_db                    = ""
params.fastani_ref                = workflow.projectDir + "/db/fastani_refs.tar.gz"
params.genome_sizes               = workflow.projectDir + "/assets/genome_sizes.json"
params.genome_references          = workflow.projectDir + "/assets/genomes.txt"

// for testing
params.sra_accessions             = []

// tool-specific command line options
params.amrfinderplus_options      = ""
params.bbduk_options              = "k=31 hdist=1"
params.bbmap_options              = ""
params.blast_db_type              = "nt"
params.blastn_options             = "-max_target_seqs 10 -max_hsps 1 -evalue 1e-25"
params.blobtools_create_options   = ""
params.blobtools_view_options     = ""
params.blobtools_plot_options     = "--format png -r species"
params.blobtools_bbmap_options    = ""
params.current_datasets           = true
params.datasets_max_genomes       = 5
params.emmtyper_options           = ''
params.extras                     = true
params.fastani_include            = true
params.fastani_options            = "--matrix"
params.fasterqdump_options        = ""
params.fastp_options              = "--detect_adapter_for_pe"
params.fastqc_options             = ""
params.iqtree2_options            = "-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000"
params.iqtree2_outgroup           = ""
params.kaptive_options            = ''
params.kleborate_options          = "-all"
params.kraken2_options            = ""
params.legsta_options             = ""
params.mash_sketch_options        = "-m 2"
params.mash_dist_options          = "-v 0 -d 0.5"
params.mash_max_hits              = 25
params.msa                        = false
params.mlst_options               = ""
params.multiqc_options            = ""
params.pbptyper_options           = ''
params.plasmidfinder_options      = ""
params.prokka_options             = "--mincontiglen 500 --compliant --locustag locus_tag --centre STAPHB"
params.quast_options              = ""
params.roary_options              = ""
params.roary_min_genes            = 1500
params.seqsero2_options           = "-m a -b mem"
params.serotypefinder_options     = ""
params.shigatyper_options         = ""
params.snp_dists_options          = "-c"
params.spades_options             = "--isolate"

// if (params.phoenix_wf) {
//   println "cp ${workflow.projectDir}/../phoenix/bin/* ${workflow.projectDir}/bin/."
//   command = ["sh", "-c", "cp ${workflow.projectDir}/../phoenix/bin/* ${workflow.projectDir}/bin/."]
//   Runtime.getRuntime().exec((String[]) command.toArray())
//   command = ["sh", "-c", "cp ${workflow.projectDir}/../phoenix/lib/* ${workflow.projectDir}/lib/."]
//   Runtime.getRuntime().exec((String[]) command.toArray())
//   include { PHOENIX_EXTERNAL } from workflow.projectDir + "/../phoenix/main" addParams(input: params.input)
//   For my future self, this doesn't work because of some of the nf-core scripts
// }

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Sharing params with subworkflows

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

include { average_nucleotide_identity } from "./subworkflows/average_nucleotide_identity"  addParams(params)
include { blobtools }                   from "./subworkflows/blobtools"                    addParams(params)
include { de_novo_alignment }           from "./subworkflows/de_novo_alignment"            addParams(params)
include { information }                 from "./subworkflows/information"                  addParams(params)
include { kraken2 }                     from "./subworkflows/kraken2"                      addParams(params)
include { min_hash_distance }           from "./subworkflows/min_hash_distance"            addParams(params)
include { phylogenetic_analysis }       from "./subworkflows/phylogenetic_analysis"        addParams(params)
include { report }                      from "./subworkflows/report"                       addParams(params)
include { test }                        from "./subworkflows/test"                         addParams(params)

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for scripts

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Creating the summary files
summary_script = Channel.fromPath(workflow.projectDir + "/bin/summary.py",     type: "file")
snpmtrx_script = Channel.fromPath(workflow.projectDir + "/bin/HeatCluster.py", type: "file")

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for input files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// using a sample sheet with the column header pf 'sample,fastq_1,fastq_2'
ch_input_reads = params.sample_sheet
  ? Channel
    .fromPath("${params.sample_sheet}", type: "file")
      .view { "Sample sheet found : ${it}" }
      .splitCsv( header: true, sep: ',' )
      .map { row -> tuple( "${row.sample}", [ file("${row.fastq_1}"), file("${row.fastq_2}") ]) }
  : Channel.empty()

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
Channel
  .fromPath("${params.gff}/*.gff", type: "file")
  .view { "gff file : $it" }
  .unique()
  .set { ch_gffs }

// Getting accession for downloading
ch_sra_accessions = Channel.from( params.sra_accessions )

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for database files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Getting the file with genome sizes of common organisms for fastqcscan. The End User can use their own file and set with a param
Channel
  .fromPath(params.genome_sizes, type: "file")
  .ifEmpty{
    println("The genome sizes file for this workflow are missing!")
    exit 1}
  .set { ch_genome_sizes }

// Getting the file with genomes to include for every run
Channel
  .fromPath(params.genome_references, type: "file")
  .ifEmpty{
    println("The genome references for this workflow are missing!")
    exit 1}
  .set { ch_genome_ref }

// Getting the reference genomes for fastANI
ch_fastani_genomes = Channel.fromPath("${params.fastani_ref}", type: "file")

// Getting the database for blobtools
ch_blast_db = params.blast_db
  ? Channel
    .fromPath(params.blast_db, type: "dir")
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
    .fromPath(params.kraken2_db, type: "dir")
    .ifEmpty{
      println("No kraken2 database was found at ${params.kraken2_db}")
      println("Set 'params.kraken2_db' to directory with kraken2 database")
      exit 1
      }
      .view { "Local kraken2 database : $it" }
  : Channel.empty()

// Getting the mash reference
ch_mash_db = params.mash_db 
  ? Channel
    .fromPath(params.mash_db, type: "file")
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
  ch_for_flag      = Channel.empty()
  ch_top_hit       = Channel.empty()

  // getting test files
  if ( ! params.sra_accessions.isEmpty() ) { 
    test(ch_sra_accessions)
    ch_raw_reads = ch_reads.mix(test.out.fastq).mix(ch_input_reads)
  } else {
    ch_raw_reads = ch_reads.mix(ch_input_reads)
  }

  // running PHOENIX first (under development)
  // if ( params.phoenix_wf )      { PHOENIX(input) }
  // running DONUT FALLS first (under development)
  // if ( params.donut_falls_wf )  { donut_falls(input) }

  // either phoenix or de_novo_alignment is required
  de_novo_alignment(ch_raw_reads)
  
  ch_contigs       = ch_fastas.mix(de_novo_alignment.out.contigs)
  ch_clean_reads   = de_novo_alignment.out.clean_reads

  ch_for_multiqc   = ch_for_multiqc.mix(de_novo_alignment.out.for_multiqc)

  // optional subworkflow blobtools (useful for interspecies contamination)
  if (params.blast_db) {
    blobtools(ch_clean_reads, ch_contigs, ch_blast_db)

    ch_for_multiqc = ch_for_multiqc.mix(blobtools.out.for_multiqc)
    ch_for_summary = ch_for_summary.mix(blobtools.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(blobtools.out.for_flag)
  }

  // optional subworkflow kraken2 (useful for interspecies contamination)
  if (params.kraken2_db) {
    kraken2(ch_clean_reads, ch_contigs, ch_kraken2_db)

    ch_for_multiqc = ch_for_multiqc.mix(kraken2.out.for_multiqc)
    ch_for_summary = ch_for_summary.mix(kraken2.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(kraken2.out.for_flag)
  } 
  
  if (params.extras) {
    // subworkflow mash for species determination
    min_hash_distance(ch_clean_reads, ch_contigs, ch_mash_db)

    ch_for_summary = ch_for_summary.mix(min_hash_distance.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(min_hash_distance.out.for_flag)

    // determining organisms in sample
    average_nucleotide_identity(
      ch_for_summary.collect(),
      ch_contigs,
      ch_fastani_genomes,
      ch_genome_ref)

    ch_for_summary = ch_for_summary.mix(average_nucleotide_identity.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(average_nucleotide_identity.out.for_flag)
    ch_top_hit     = ch_top_hit.mix(average_nucleotide_identity.out.top_hit)
    ch_datasets    = average_nucleotide_identity.out.datasets_summary.ifEmpty('none')

    ch_top_hit_files = ch_top_hit.map {it -> tuple(it[0], it[1])}
    ch_top_hit_hit   = ch_top_hit.map {it -> tuple(it[0], it[2])}

    ch_contigs
      .join(min_hash_distance.out.mash_err)
      .join(min_hash_distance.out.for_flag)
      .join(average_nucleotide_identity.out.for_flag)
      .combine(ch_genome_sizes)
      .combine(ch_datasets)
      .set{ ch_size }

    // getting all the other information
    information(
      ch_raw_reads, 
      ch_contigs, 
      ch_for_flag, 
      ch_size)

    ch_for_summary = ch_for_summary.mix(information.out.for_summary)
    ch_for_multiqc = ch_for_multiqc.mix(information.out.for_multiqc)
  } 

  // optional subworkflow for comparing shared genes
  if ( params.msa ) { 
    phylogenetic_analysis(
      snpmtrx_script,
      ch_contigs.ifEmpty([]), 
      ch_gffs.ifEmpty([]), 
      ch_top_hit.ifEmpty([]))
    
    ch_for_multiqc   = ch_for_multiqc.mix(phylogenetic_analysis.out.for_multiqc)
  }

  // getting a summary of everything
  if (params.extras) { 
    report(
      ch_raw_reads, 
      ch_fastas, 
      ch_for_multiqc.collect(), 
      ch_for_summary.concat(summary_script).collect())
  }
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
