#!/usr/bin/env nextflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Welcome to this workflow! Issues and contributions are gladly accepted at https://github.com/UPHL-BioNGS/Grandeur .

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

println("Currently using the Grandeur workflow for use with microbial sequencing. The view is great from 8299 feet (2530 meters) above sea level.\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: ${workflow.manifest.version}")
println("")

nextflow.enable.dsl    = 2
nextflow.enable.strict = true

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Getting config file

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

params.config_file  = false
if ( params.config_file ) {
  def src = new File("${workflow.projectDir}/configs/grandeur_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println("A config file can be found at ${workflow.launchDir}/edit_me.config")
  exit 0
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Defining params

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

params.outdir               = "grandeur"

// input files
params.reads                = ""
params.fastas               = ""
params.sample_sheet         = ""
params.fasta_list           = ""

// external files
params.kraken2_db           = ""
params.blast_db             = ""
params.mash_db              = ""
params.fastani_ref          = ""
params.fastani_ref_list     = ""
params.genome_sizes         = workflow.projectDir + "/assets/genome_sizes.json"

// for downloading from databases
params.sra_accessions       = []

// thresholds and other params
params.minimum_reads        = 10000
params.datasets_max_genomes = 5
params.mash_max_hits        = 25
params.min_core_genes       = 1500
params.iqtree2_outgroup     = ""

// subworkflow flags
params.current_datasets     = false
params.skip_extras          = false
params.exclude_top_hit      = false
params.msa                  = false


// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Checking params

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

def paramCheck(keys) {
  set_keys = [
    "outdir",
    "fastas",
    "msa",
    "kraken2_db",
    "mash_db",
    "config_file",
    "reads",
    "sample_sheet",
    "fasta_list",
    "blast_db",
    "fastani_ref",
    "fastani_ref_list",
    "iqtree2_outgroup",
    "genome_sizes",
    "sra_accessions",
    "minimum_reads",
    "datasets_max_genomes",
    "mash_max_hits",
    "min_core_genes",
    "current_datasets",
    "skip_extras",
    "exclude_top_hit"]

  for(key in keys){
    if (key !in set_keys){
      println("FATAL: ${key} isn't a supported param!")
      println("Supported params: ${set_keys}")
      exit 1
    }
  }
}

paramCheck(params.keySet())

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Sharing params with subworkflows

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// TODO : https://oldsite.nf-co.re/pipeline_schema_builder nf-core schema build

include { average_nucleotide_identity }   from "./subworkflows/average_nucleotide_identity"   addParams(params)
//include { blobtools }                     from "./subworkflows/blobtools"                     addParams(params)
include { de_novo_alignment }             from "./subworkflows/de_novo_alignment"             addParams(params)
include { information }                   from "./subworkflows/information"                   addParams(params)
include { kmer_taxonomic_classification } from "./subworkflows/kmer_taxonomic_classification" addParams(params)
include { min_hash }                      from "./subworkflows/min_hash"                      addParams(params)
include { phylogenetic_analysis }         from "./subworkflows/phylogenetic_analysis"         addParams(params)
include { report }                        from "./subworkflows/report"                        addParams(params)
include { test }                          from "./subworkflows/test"                          addParams(params)

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for scripts

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

dataset_script = Channel.fromPath(workflow.projectDir + "/bin/datasets_download.py", type: "file")
evaluat_script = Channel.fromPath(workflow.projectDir + "/bin/evaluate.py",          type: "file")
jsoncon_script = Channel.fromPath(workflow.projectDir + "/bin/json_convert.py",      type: "file")
multiqc_script = Channel.fromPath(workflow.projectDir + "/bin/for_multiqc.py",       type: "file")
summary_script = Channel.fromPath(workflow.projectDir + "/bin/summary.py",           type: "file")
summfle_script = Channel.fromPath(workflow.projectDir + "/bin/summary_file.py",      type: "file")
version_script = Channel.fromPath(workflow.projectDir + "/bin/versions.py",          type: "file")

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for input files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

if (params.sample_sheet) {
  // using a sample sheet with the column header of 'sample,fastq_1,fastq_2'
  Channel
    .fromPath("${params.sample_sheet}", type: "file")
    .view { "Sample sheet found : ${it}" }
    .splitCsv( header: true, sep: ',' )
    .map { row ->
      meta = [id:row.sample]
      tuple( meta, [ 
        file("${row.fastq_1}", checkIfExists: true), 
        file("${row.fastq_2}", checkIfExists: true)])
    }
    .set {ch_reads}

} else {
  // Getting the fastq files from a directory
  ch_reads = params.reads
    ? Channel
        .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                        "${params.reads}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
        .map { it ->
          meta = [id:it[0].replaceAll(~/_S[0-9]+_L[0-9]+/,"")] 
          tuple( meta, 
            file(it[0], checkIfExists: true), 
            file(it[1], checkIfExists: true))
        }
        .unique()
        .view { "Paired-end fastq files found : ${it[0]}" }
    : Channel.empty()
}

if (params.fasta_list) {
  // TODO : make sure this works

  // getting fastas from a file
  Channel
    .fromPath("${params.fasta_list}", type: "file")
    .view { "Fasta list found : ${it}" }
    .splitText()
    .map{ it -> it.trim()}
    .map{ it -> file(it) }
    .map { it ->
      meta = [id:it.baseName]
      tuple( meta, it)
    }
    .set{ ch_fastas }
} else {
  // getting fastas from a directory
  ch_fastas = params.fastas
    ? Channel
      .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
      .view { "fasta file : $it" }
      .map { it ->
        meta = [id: it.baseName]
        tuple( meta, file(it, checkIfExists: true))
      }
      .unique()
    : Channel.empty()
}

// Getting accession for downloading

// from SRA
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

//# user supplied fastani reference genomes
ch_fastani_genomes = Channel.empty()

if ( params.fastani_ref ) {

  // TODO : test fastani_ref
  Channel
    .of( params.fastani_ref )
    .splitCsv()
    .flatten()
    // no meta id
    .map { it -> file(it) }
    .view{ "Additional fastani reference genomes : $it" }
    .set { ch_fastani_genomes_input }

  ch_fastani_genomes = ch_fastani_genomes.mix(ch_fastani_genomes_input)
}

if ( params.fastani_ref_list ) {
  Channel.fromPath(params.fastani_ref_list, type: "file")
    .splitText()
    .map( it -> it.trim())
    .map{ it -> file(it) }
    .view{ "Additional fastani reference genome from file : $it" }
    .set{ ch_fastani_ref_list }

  ch_fastani_genomes = ch_fastani_genomes.mix(ch_fastani_ref_list)
}

println("The files and directory for results is " + params.outdir )

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Workflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow {
  ch_for_multiqc  = Channel.empty()
  ch_for_summary  = ch_genome_sizes
  ch_for_flag     = Channel.empty()
  ch_versions     = Channel.empty()

  // getting test files
  if ( ! params.sra_accessions.isEmpty() ) { 
    test(ch_sra_accessions)
    ch_raw_reads = ch_reads.mix(test.out.fastq)
    ch_versions = ch_versions.mix(test.out.versions)
  } else {
    ch_raw_reads = ch_reads
  }

  if ( params.sample_sheet || params.reads ) {
    de_novo_alignment(ch_raw_reads)

    ch_contigs     = ch_fastas.mix(de_novo_alignment.out.contigs)
    ch_clean_reads = de_novo_alignment.out.clean_reads
    ch_for_multiqc = ch_for_multiqc.mix(de_novo_alignment.out.for_multiqc)
    ch_versions = ch_versions.mix(de_novo_alignment.out.versions)

  } else {
    ch_contigs     = ch_fastas
    ch_clean_reads = Channel.empty()
  }

  // optional subworkflow blobtools (useful for interspecies contamination)
  // if ( params.blast_db && ( params.sample_sheet || params.reads )) {
  //   blobtools(ch_clean_reads, ch_contigs, ch_blast_db )

  //   ch_for_multiqc = ch_for_multiqc.mix(blobtools.out.for_multiqc)
  //   ch_for_summary = ch_for_summary.mix(blobtools.out.for_summary)
  //   ch_for_flag    = ch_for_flag.mix(blobtools.out.for_flag)
  // ch_versions = ch_versions.mix(blobtools.out.versions)
  // }

  // optional subworkflow kraken2 (useful for interspecies contamination)
  if ( params.kraken2_db && ( params.sample_sheet || params.reads )) {
    // TODO : figure out where the kraken2 fastq files are

    kmer_taxonomic_classification(ch_clean_reads, ch_kraken2_db )

    ch_for_multiqc = ch_for_multiqc.mix(kmer_taxonomic_classification.out.for_multiqc)
    ch_for_summary = ch_for_summary.mix(kmer_taxonomic_classification.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(kmer_taxonomic_classification.out.for_flag)
    ch_versions    = ch_versions.mix(kmer_taxonomic_classification.out.versions)
  } 
  
  if ( ! params.skip_extras ) {
    // subworkflow mash for species determination
    min_hash(ch_clean_reads, ch_fastas, ch_mash_db)

    // determining organisms in sample
    average_nucleotide_identity(
      ch_for_summary.mix(min_hash.out.for_summary).collect(),
      ch_contigs,
      ch_fastani_genomes.ifEmpty([]),
      dataset_script)

    ch_for_flag = ch_for_flag.mix(average_nucleotide_identity.out.for_flag).mix(min_hash.out.for_flag)
    ch_top_hit  = average_nucleotide_identity.out.top_hit

    // getting all the other information
    information(
      ch_raw_reads, 
      ch_contigs, 
      ch_for_flag, 
      summfle_script,
      jsoncon_script)

    ch_for_summary = ch_for_summary.mix(information.out.for_summary).mix(min_hash.out.for_summary).mix(average_nucleotide_identity.out.for_summary)
    ch_for_multiqc = ch_for_multiqc.mix(information.out.for_multiqc)
    ch_versions    = ch_versions.mix(min_hash.out.versions).mix(average_nucleotide_identity.out.versions).mix(information.out.versions)
  } else {
    ch_top_hit = Channel.empty()
  }

  // optional subworkflow for comparing shared genes
  if ( params.msa ) {
    phylogenetic_analysis(
      evaluat_script,
      ch_contigs.ifEmpty([]),
      ch_top_hit.ifEmpty([]))
    
    ch_for_multiqc = ch_for_multiqc.mix(phylogenetic_analysis.out.for_multiqc)
    ch_versions    = ch_versions.mix(phylogenetic_analysis.out.versions)
  }

  // getting a summary of everything
  if ( ! params.skip_extras ) {
    report(
      ch_raw_reads, 
      ch_fastas, 
      ch_for_multiqc.collect(), 
      ch_for_summary.concat(summary_script).collect(),
      ch_versions.collect(),
      multiqc_script,
      version_script)
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
