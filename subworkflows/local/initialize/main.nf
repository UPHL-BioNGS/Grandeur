include { TEST } from "../../subworkflows/local/test"

def paramCheck(keys) {
  def set_keys = [
    "outdir",
    "input",
    "fastas",
    "msa",
    "kraken2_db",
    "mash_db",
    "config_file",
    "reads",
    "sample_sheet",
    "fasta_list",
    "blast_db",
    "blast_db_type",
    "fastani_ref",
    "fastani_ref_list",
    "iqtree2_outgroup",
    "genome_sizes",
    "sra_accessions",
    "genome_accessions",
    "minimum_reads",
    "datasets_max_genomes",
    "mash_max_hits",
    "min_core_genes",
    "min_core_per",
    "current_datasets",
    "annotator",
    "skip_extras",
    "exclude_top_hit",
    "aligner",
    "publish_dir_mode",
    "email",
    "email_on_fail",
    "plaintext_email",
    "monochrome_logs",
    "hook_url",
    "help",
    "version",
    "pipelines_testdata_base_path",
    "config_profile_name",
    "config_profile_description",
    "custom_config_version",
    "custom_config_base",
    "config_profile_contact",
    "config_profile_url",
    "validation-fail-unrecognised-params",
    "validationFailUnrecognisedParams",
    "validation-lenient-mode",
    "validationLenientMode",
    "validationShowHiddenParams",
    "validation-show-hidden-params",
    "validate_params"
  ]

  keys.each { x ->
    if (x !in set_keys){
      println("WARNING: ${x} isn't a supported param!")
      println("Supported params: ${set_keys}")
    }
  }
}

workflow INITIALIZE {
  main:
  ch_fastas    = Channel.empty()
  ch_versions  = Channel.empty()
  
  //# For aesthetics - and, yes, we are aware that there are better ways to write this than a bunch of 'println' statements
  println('') 
  println('   /^^^^    /^^^^^^^           /^        /^^^     /^^ /^^^^^     /^^^^^^^^ /^^     /^^ /^^^^^^^    ')
  println(' /^    /^^  /^^    /^^        /^ ^^      /^ /^^   /^^ /^^   /^^  /^^       /^^     /^^ /^^    /^^  ')
  println('/^^         /^^    /^^       /^  /^^     /^^ /^^  /^^ /^^    /^^ /^^       /^^     /^^ /^^    /^^  ')
  println('/^^         /^ /^^          /^^   /^^    /^^  /^^ /^^ /^^    /^^ /^^^^^^   /^^     /^^ /^ /^^      ')
  println('/^^   /^^^^ /^^  /^^       /^^^^^^ /^^   /^^   /^ /^^ /^^    /^^ /^^       /^^     /^^ /^^  /^^    ')
  println(' /^^    /^  /^^    /^^    /^^       /^^  /^^    /^ ^^ /^^   /^^  /^^       /^^     /^^ /^^    /^^  ')
  println('  /^^^^^    /^^      /^^ /^^         /^^ /^^      /^^ /^^^^^     /^^^^^^^^   /^^^^^    /^^      /^^')
  println('')                                                                            

  println("Currently using the Grandeur workflow for use with microbial sequencing.")
  println("The view is great from 8299 feet (2530 meters) above sea level.\n")
  println("Author: Erin Young")
  println("email: eriny@utah.gov")
  println("Version: ${workflow.manifest.version}")
  println("")


  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  // Getting config file

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  params.config_file  = false
  if ( params.config_file ) {
    def src = new File("${workflow.projectDir}/configs/grandeur_template.config")
    def dst = new File("${workflow.launchDir}/edit_me.config")
    dst << src.text
    println("A config file can be found at ${workflow.launchDir}/edit_me.config")

    def src1 = new File("${workflow.projectDir}/configs/grandeur_params.yml")
    def dst1 = new File("${workflow.launchDir}/edit_me.yml")
    dst1 << src1.text
    println("A params file can be found at ${workflow.launchDir}/edit_me.yml")
    exit 0
  }

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  // Checking params

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  paramCheck(params.keySet())



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
        def meta = [id:row.sample]
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
            def meta = [id:it[0].replaceAll(~/_S[0-9]+_L[0-9]+/,"")] 
            tuple( meta, [
              file(it[1][0], checkIfExists: true), 
              file(it[1][1], checkIfExists: true)])
          }
          .unique()
          .view { "Paired-end fastq files found : ${it[0].id}" }
      : Channel.empty()
  }

  if (params.fasta_list) {
    // getting fastas from a file
    Channel
      .fromPath("${params.fasta_list}", type: "file")
      .view { "Fasta list found : ${it}" }
      .splitText()
      .map{ it -> it.trim()}
      .map{ it -> file(it) }
      .map { it ->
        def meta = [id:it.baseName]
        tuple( meta, it)
      }
      .set{ ch_fastas }
  } else {
    // getting fastas from a directory
    ch_fastas = params.fastas
      ? Channel
        .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
        .view { "Fasta file found : ${it.baseName}" }
        .map { it ->
          def meta = [id: it.baseName]
          tuple( meta, file(it, checkIfExists: true))
        }
        .unique()
      : Channel.empty()
  }

  // Getting accession for downloading

  // from SRA
  ch_sra_accessions   = Channel.from( params.sra_accessions )

  // from genomes
  ch_genome_accessions = Channel.from( params.genome_accessions)

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
      .map{ it -> it.trim()}
      .map{ it -> file(it) }
      .view{ "Additional fastani reference genome from file : $it" }
      .set{ ch_fastani_ref_list }

    ch_fastani_genomes = ch_fastani_genomes.mix(ch_fastani_ref_list)
  }

  println("The files and directory for results is " + params.outdir )

  // getting test files
  if ( ! params.sra_accessions.isEmpty()  || ! params.genome_accessions.isEmpty() ) { 
    TEST(
      ch_sra_accessions.ifEmpty([]), 
      ch_genome_accessions.ifEmpty([])
    )
    ch_reads    = ch_reads.mix(TEST.out.fastq)
    ch_fastas   = ch_fastas.mix(TEST.out.fasta)
    ch_versions = TEST.out.versions
  }

  emit:
  reads           = ch_reads
  fastas          = ch_fastas
  fastani_genomes = ch_fastani_genomes
  versions        = ch_versions
  genome_sizes    = ch_genome_sizes
  mash_db         = ch_mash_db
  kraken2_db      = ch_kraken2_db
  blast_db        = ch_blast_db
  dataset_script  = dataset_script
  evaluat_script  = evaluat_script
  jsoncon_script  = jsoncon_script
  multiqc_script  = multiqc_script
  summary_script  = summary_script
  summfle_script  = summfle_script
  version_script  = version_script
}