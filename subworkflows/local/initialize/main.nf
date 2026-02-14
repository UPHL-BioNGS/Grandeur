include { TEST } from "../../../subworkflows/local/test"

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
    "checkm2_db",
    "sylph_db",
    "reference_genomes",
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
  ch_fastas    = channel.empty()
  ch_versions  = channel.empty()
  
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

  // channels for scripts

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  dataset_script = channel.fromPath(workflow.projectDir + "/bin/datasets_download.py", type: "file")
  evaluat_script = channel.fromPath(workflow.projectDir + "/bin/evaluate.py",          type: "file")
  jsoncon_script = channel.fromPath(workflow.projectDir + "/bin/json_convert.py",      type: "file")
  multiqc_script = channel.fromPath(workflow.projectDir + "/bin/for_multiqc.py",       type: "file")
  summary_script = channel.fromPath(workflow.projectDir + "/bin/summary.py",           type: "file")
  summfle_script = channel.fromPath(workflow.projectDir + "/bin/summary_file.py",      type: "file")
  version_script = channel.fromPath(workflow.projectDir + "/bin/versions.py",          type: "file")

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  // channels for input files

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  if (params.sample_sheet) {
    // using a sample sheet with the column header of 'sample,fastq_1,fastq_2'
    channel
      .fromPath("${params.sample_sheet}", type: "file")
      .view { "Sample sheet found : ${it}" }
      .splitCsv( header: true, sep: ',' )
      .map { row ->
        def meta = [id:row.sample]
        tuple( meta, [ 
          file("${row.fastq_1}", checkIfExists: true), 
          file("${row.fastq_2}", checkIfExists: true)])
      }
      .unique()
      .set {ch_reads}

  } else {
    // Getting the fastq files from a directory
    ch_reads = params.reads
      ? channel
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
      : channel.empty()
  }

  if (params.fasta_list) {
    // getting fastas from a file
    channel
      .fromPath("${params.fasta_list}", type: "file", checkIfExists : true)
      .view { "Fasta list found : ${it}" }
      .splitText()
      .map{ it -> it.trim()}
      .map{ it -> file(it) }
      .map { it ->
        def meta = [id:it.baseName]
        tuple( meta, it)
      }
      .unique()
      .set{ ch_fastas }
  } else {
    // getting fastas from a directory
    ch_fastas = params.fastas
      ? channel
        .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
        .view { "Fasta file found : ${it.baseName}" }
        .map { it ->
          def meta = [id: it.baseName]
          tuple( meta, file(it, checkIfExists: true))
        }
        .unique()
      : channel.empty()
  }

  // Getting accession for downloading

  // from SRA
  ch_sra_accessions   = channel.from( params.sra_accessions )

  // from genomes
  ch_genome_accessions = channel.from( params.genome_accessions)

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  // channels for database files

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  // Getting the file with genome sizes of common organisms for fastqcscan. The End User can use their own file and set with a param
  if ( params.genome_sizes ) {
    println("Using genome sizes file at ${params.genome_sizes}")
    channel
      .fromPath(params.genome_sizes, type: "file", checkIfExists: true)
      .ifEmpty{
        println("The custom genome sizes file for this workflow are missing!")
        exit 1}
      .set { ch_genome_sizes }
  } else {
    channel
      .fromPath("${projectDir}/assets/genome_sizes.json", type: "file", checkIfExists: true)
      .ifEmpty{
        println("The genome sizes file for this workflow are missing!")
        exit 1}
      .set { ch_genome_sizes }
  }

  // Getting the kraken2 database
  ch_kraken2_db = params.kraken2_db
    ? channel
      .fromPath(params.kraken2_db, type: "dir")
      .ifEmpty{
        println("No kraken2 database was found at ${params.kraken2_db}")
        println("Set 'params.kraken2_db' to directory with kraken2 database")
        exit 1
        }
        .view { "Local kraken2 database : $it" }
    : channel.empty()

  // Getting the mash reference
  ch_mash_db = params.mash_db 
    ? channel
      .fromPath(params.mash_db, type: "file")
      .ifEmpty{
        println("No mash database was found at ${params.mash_db}")
        println("Set 'params.mash_db' to file of pre-sketched mash reference")
        exit 1
        }
      .view { "Mash reference : $it" }
    : channel.empty()

  // Getting the kraken2 database
  ch_checkm2_db = params.checkm2_db
    ? channel
      .fromPath(params.checkm2_db, type: "dir")
      .ifEmpty{
        println("No checkm2 database was found at ${params.checkm2_db}")
        println("Set 'params.checkm2_db' to directory with checkm2 database")
        exit 1
        }
        .view { "Local checkm2 database : $it" }
    : channel.empty()


  ch_sylph_db = params.sylph_db
    ? channel
      .fromPath(params.sylph_db, type: "dir")
      .ifEmpty{
        println("No Sylph database was found at ${params.sylph_db}")
        println("Set 'params.sylph_db' to directory with Sylph database")
        exit 1
        }
        .view { "Local Sylph database : $it" }
    : channel.empty()
  
  // if using additional fasta files for skani
  if (  params.reference_genomes ) {
    channel.fromPath(params.reference_genomes, type: "file")
      .splitText()
      .map{ it -> it.trim()}
      .map{ it -> file(it) }
      .unique()
      .view{ "Additional reference genome from file : $it" }
      .set{ ch_reference_genomes }
  } else {
    ch_reference_genomes = channel.empty()
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
  reads             = ch_reads
  fastas            = ch_fastas
  reference_genomes = ch_reference_genomes
  versions          = ch_versions
  genome_sizes      = ch_genome_sizes
  mash_db           = ch_mash_db
  kraken2_db        = ch_kraken2_db
  checkm2_db        = ch_checkm2_db
  sylph_db          = ch_sylph_db
  dataset_script    = dataset_script
  evaluat_script    = evaluat_script
  jsoncon_script    = jsoncon_script
  multiqc_script    = multiqc_script
  summary_script    = summary_script
  summfle_script    = summfle_script
  version_script    = version_script


}

workflow.onComplete {
  println("Inititalization completed at: $workflow.complete")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}