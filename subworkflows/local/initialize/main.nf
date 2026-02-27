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
      log.warn "WARNING: ${x} isn't a supported param!"
      log.info "Supported params: ${set_keys}"
    }
  }
}

workflow INITIALIZE {
  main:
  ch_versions  = channel.empty()
  
  log.info """\

   /^^^^    /^^^^^^^           /^        /^^^     /^^ /^^^^^     /^^^^^^^^ /^^     /^^ /^^^^^^^    
 /^    /^^  /^^    /^^        /^ ^^      /^ /^^   /^^ /^^   /^^  /^^       /^^     /^^ /^^    /^^  
/^^         /^^    /^^       /^  /^^     /^^ /^^  /^^ /^^    /^^ /^^       /^^     /^^ /^^    /^^  
/^^         /^ /^^          /^^   /^^    /^^  /^^ /^^ /^^    /^^ /^^^^^^   /^^     /^^ /^ /^^      
/^^   /^^^^ /^^  /^^       /^^^^^^ /^^   /^^   /^ /^^ /^^    /^^ /^^       /^^     /^^ /^^  /^^    
 /^^    /^  /^^    /^^    /^^       /^^  /^^    /^ ^^ /^^   /^^  /^^       /^^     /^^ /^^    /^^  
  /^^^^^    /^^      /^^ /^^         /^^ /^^      /^^ /^^^^^     /^^^^^^^^   /^^^^^    /^^      /^^

Currently using the Grandeur workflow for use with microbial sequencing.
The view is great from 8299 feet (2530 meters) above sea level.

Author: Erin Young
email: eriny@utah.gov
Version: ${workflow.manifest.version}
"""

log.info """
------------------------------------------------------------------------------------------------------------

Initializing Workflow and Evaluating Parameters

------------------------------------------------------

"""

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  // Getting config file

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  params.config_file  = false
  if ( params.config_file ) {
    def src = new File("${workflow.projectDir}/configs/grandeur_template.config")
    def dst = new File("${workflow.launchDir}/edit_me.config")
    dst << src.text
    log.info "A config file can be found at ${workflow.launchDir}/edit_me.config"

    def src1 = new File("${workflow.projectDir}/configs/grandeur_params.yml")
    def dst1 = new File("${workflow.launchDir}/edit_me.yml")
    dst1 << src1.text
    log.info "A params file can be found at ${workflow.launchDir}/edit_me.yml"
    exit 0
  }

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  // Checking params

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  paramCheck(params.keySet())

  log.info "Documentation for this workflow can be found at https://github.com/UPHL-BioNGS/Grandeur/wiki\n"
  log.info "All files will be saved to ${params.outdir}\n\t- To change this, set 'params.outdir' to the desired output directory.\n"


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

  log.info """
------------------------------------------------------

Initializing Sample Input Files

------------------------------------------------------

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ param             ┃ type  ┃ value                              
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ sample_sheet      ┃ file  ┃ ${params.sample_sheet}
┃ reads             ┃ dir   ┃ ${params.reads}
┃ fasta_list        ┃ file  ┃ ${params.fasta_list}
┃ fastas            ┃ dir   ┃ ${params.fastas}
┃ sra_accessions    ┃ list  ┃ ${params.sra_accessions}
┃ genome_accessions ┃ list  ┃ ${params.sra_accessions}
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━

------------------------------------------------------
"""


  if (params.sample_sheet) {
    log.info "Using sample sheet at ${params.sample_sheet}"
    log.info "\t- The sample sheet should be a csv file with the column header of 'sample,fastq_1,fastq_2' and the respective values for each sample listed below. The base name of each sample (the value in the 'sample' column) is used as the \"meta.id\" value, and is used when generating output files and summarizing results."
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
      .ifEmpty{
        log.fatal "The 'params.sample_sheet' was set, but no input files were found!"
        exit 1}
      .set {ch_reads}

  } else {
    // Getting the FASTQ files from a directory
    if (params.reads) {
      log.info "Looking for FASTQ files in directory ${params.reads}"
      log.info "\t- FASTQ files should have the extension .fastq, .fastq.gz, .fq, or .fq.gz"
      channel
          .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                          "${params.reads}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
          .map { it ->
            def meta = [id:it[0].replaceAll(~/_S[0-9]+_L[0-9]+/,"")] 
            tuple( meta, [
              file(it[1][0], checkIfExists: true), 
              file(it[1][1], checkIfExists: true)])
          }
          .unique()
          .view { "Paired-end FASTQ files found : ${it[0].id}" }
          .ifEmpty{
            log.fatal "The 'params.reads' was set, but no input files were found!"
            exit 1}
          .set { ch_reads }
    } else {
      log.info "FYI: Input FASTQ files can be provided to Grandeur with 'params.reads' or with a sample sheet designated with 'params.sample_sheet'."
      ch_reads = channel.empty()
    }
  }

  if (params.fasta_list) {
    log.info "Loading FASTA list at ${params.fasta_list}"
    // getting FASTAs from a file
    channel
      .fromPath("${params.fasta_list}", type: "file", checkIfExists : true)
      .view { "FASTA list found : ${it}" }
      .splitText()
      .map{ it -> it.trim()}
      .map{ it -> file(it) }
      .map { it ->
        def meta = [id:it.baseName]
        tuple( meta, it)
      }
      .unique()
      .ifEmpty{
          log.fatal "The 'params.fasta_list' was set, but no input files were found!"
          exit 1}
      .set{ ch_fastas }
  } else {
    // getting FASTAs from a directory
    if (params.fastas) {
      log.info "Looking for FASTA files in directory ${params.fastas}"
      log.info "\t- FASTA files should have the extension .fa, .fasta, or .fna"
      log.info "\t- The base name of each FASTA file is used as the \"meta.id\" value, and is used when generating output files and summarizing results."
      channel
        .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
        .view { "FASTA file found : ${it.baseName}" }
        .map { it ->
          def meta = [id: it.baseName]
          tuple( meta, file(it, checkIfExists: true))
        }
        .unique()
        .ifEmpty{
          log.fatal "The 'params.fastas' was set, but no input files were found!"
          exit 1}
        .set { ch_fastas }
    } else {
      log.info "FYI: Input FASTA files can be provided to Grandeur with 'params.fastas' or with a list of FASTA files designated with 'params.fasta_list'."
      ch_fastas = channel.empty()
    }
  }

  // Getting accession for downloading

  // from SRA
  if (params.sra_accessions) {
    log.info "Loading SRA accessions listed in ${params.sra_accessions}"
    channel
      .from(params.sra_accessions)
      .filter{ it -> it }
      .unique()
      .view { "Using SRA accession : ${it}" }
      .ifEmpty{
        log.fatal "The 'params.sra_accessions' was set, but no value was given!"
        exit 1}
      .set { ch_sra_accessions }
  } else {
    ch_sra_accessions = channel.empty()
  }

  // from genomes
  if (params.genome_accessions) {
    log.info "Loading genome accessions listed in ${params.genome_accessions}"
    channel
      .from(params.genome_accessions)
      .filter{ it -> it }
      .unique()
      .view { "Using Genome accession : ${it}" }
      .ifEmpty{
        log.fatal "The 'params.genome_accessions' was set, but no value was given!"
        exit 1}
      .set { ch_genome_accessions }
  } else {
    ch_genome_accessions = channel.empty()
  }

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  // channels for database files

  // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

  log.info """
------------------------------------------------------

Initializing Databases and References

------------------------------------------------------

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ param             ┃ type  ┃ value                              
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ kraken2_db        ┃ dir   ┃ ${params.kraken2_db}
┃ mash_db           ┃ file  ┃ ${params.mash_db}
┃ checkm2_db        ┃ file  ┃ ${params.checkm2_db}
┃ sylph_db          ┃ file  ┃ ${params.sylph_db}
┃ reference_genomes ┃ file  ┃ ${params.reference_genomes}
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━

"""


  // Getting the file with genome sizes of common organisms for summary.
  if ( params.genome_sizes ) {
    log.info "Using custom genome sizes file at ${params.genome_sizes}"
    log.warn "\tThis is most-often used for dev and testing purposes."
    channel
      .fromPath(params.genome_sizes, type: "file", checkIfExists: true)
      .ifEmpty{
        log.info "The custom genome sizes file for this workflow are missing!"
        exit 1}
      .set { ch_genome_sizes }
  } else {
    channel
      .fromPath("${projectDir}/assets/genome_sizes.json", type: "file", checkIfExists: true)
      .ifEmpty{
        log.info "The genome sizes file for this workflow are missing!"
        exit 1}
      .set { ch_genome_sizes }
  }

  // Getting the kraken2 database
  if (params.kraken2_db){
    log.info "Looking for KRAKEN2 database in directory ${params.kraken2_db}"
    channel
      .fromPath(params.kraken2_db, type: "dir")
      .ifEmpty{
        log.info "No KRAKEN2 database was found at ${params.kraken2_db}"
        log.info "Set 'params.kraken2_db' to **directory** with KRAKEN2 database"
        exit 1
        }
        .view { "Using KRAKEN2 database : $it" }
        .set { ch_kraken2_db }
  } else {
    log.info "FYI: A KRAKEN2 database can be loaded into Grandeur with 'params.kraken2_db'."
    ch_kraken2_db = channel.empty()
  }

  // Getting the mash reference
  if (params.mash_db) {
    log.info "Looking for MASH database file at ${params.mash_db}"
    channel
      .fromPath(params.mash_db, type: "file", checkIfExists: true)
      .ifEmpty{
        log.info "No MASH database was found at ${params.mash_db}"
        log.info "Set 'params.mash_db' to file of pre-sketched MASH reference"
        exit 1
        }
        .view { "Using MASH reference : $it" }
        .set { ch_mash_db }
  } else {
    log.info "Using default MASH database located in STaPH-B/mash container (RefSeqSketchesDefaults.msh)."
    log.info "FYI: A custome MASH database can be loaded into Grandeur with 'params.mash_db', more information can be found at https://github.com/UPHL-BioNGS/Grandeur/wiki/mash"
    ch_mash_db = channel.empty()
  }

  // Getting the kraken2 database
  if (params.checkm2_db){
    log.info "Looking for CHECKM2 database file at ${params.checkm2_db}"
    channel
      .fromPath(params.checkm2_db, type: "dir")
      .ifEmpty{
        log.info "No CHECKM2 database was found at ${params.checkm2_db}"
        log.info "Set 'params.checkm2_db' to **directory** with CHECKM2 database"
        exit 1
        }
        .view { "Using CHECKM2 database : $it" }
        .set { ch_checkm2_db }
  } else {
    log.info "FYI: A CHECKM2 database can be loaded into Grandeur with 'params.checkm2_db'."
    log.info "\t- Please read the wiki for instructions on how to create a CHECKM2 database for use with Grandeur at https://github.com/UPHL-BioNGS/Grandeur/wiki/checkm2_db"
    ch_checkm2_db = channel.empty()
  }

  // Getting the sylph database
  if (params.sylph_db){
    log.info "Looking for SYLPH database in directory ${params.sylph_db}"
    channel
      .fromPath(params.sylph_db, type: "dir")
      .ifEmpty{
        log.info "No SYLPH database was found at ${params.sylph_db}"
        log.info "Set 'params.sylph_db' to **directory** with SYLPH database"
        exit 1
        }
        .view { "Using SYLPH database : $it" }
        .set { ch_sylph_db }
  } else {
    log.info "FYI: A SYLPH database can be loaded into Grandeur with 'params.sylph_db'."
    log.info "\t- Please read the wiki for instructions on how to create a SYLPH database for use with Grandeur at https://github.com/UPHL-BioNGS/Grandeur/wiki/sylph_db"
    ch_sylph_db = channel.empty()
  }
  
  // if using additional fasta files for ani
  if (  params.reference_genomes ) {
    log.info "Loading additional reference genomes listed in ${params.reference_genomes}"
    log.info "\t- Please note that these files must be named genus_species_uniquename.fasta."
    channel.fromPath(params.reference_genomes, type: "file")
      .splitText()
      .map{ it -> it.trim()}
      .map{ it -> file(it) }
      .unique()
      .view{ "Additional reference genome from file : $it" }
      .set{ ch_reference_genomes }
  } else {
    log.info "FYI: Additional reference genomes can be loaded into Grandeur for ANI analysis with 'params.reference_genomes'."
    log.info "\t- Please note that this file should list the path for one reference genome per line, and these references must be named genus_species_uniquename.fasta."
    ch_reference_genomes = channel.empty()
  }

  log.info """
------------------------------------------------------

Initializing Workflow Options

------------------------------------------------------

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ param             ┃ type  ┃ value                              
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ msa               ┃ bool  ┃ ${params.msa}
┃ skip_extras       ┃ bool  ┃ ${params.skip_extras}
┃ aligner           ┃ str   ┃ ${params.aligner}
┃ minimum_reads     ┃ int   ┃ ${params.minimum_reads}
┃ min_core_genes    ┃ int   ┃ ${params.min_core_genes}
┃ min_core_per      ┃ float ┃ ${params.min_core_per}
┃ current_datasets  ┃ bool  ┃ ${params.current_datasets}
┃ exclude_top_hit   ┃ bool  ┃ ${params.exclude_top_hit}
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━


"""


  if (params.msa) {
    log.info "'params.msa' is set to true. All input files will be put through the PHYLOGENETIC_ANALYSIS subworkflow. Please ensure that all input genomes are reasonably related."
    log.info "The minimum number of genes these genomes should share is ${params.min_core_genes} and the minimum core genome percentage is ${params.min_core_per}%. These can be adjusted with 'params.min_core_genes' and 'params.min_core_per'."
  } else {
    log.info "FYI: The PHYLOGENETIC_ANALYSIS subworkflow is skipped by default. To compare isolates with reasonable top hits, set 'params.msa' to true."
  }


  if ( params.skip_extras ) {
    log.info "'params.skip_extras' is set to true. Skipping all \"extra\" processes and subworkflows. This focuses on the core assembly of reads (if FASTQ files are provided) or multiple sequence alignment (if 'params.msa' is set to true)."
  } else {
    log.info "FYI: It is possible to skip ANI analysis, subtyping, and taxonomic profiling subworkflows. To skip these steps, set 'params.skip_extras' to true."
  }

  // getting test files
  if ( ! params.sra_accessions.isEmpty()  || ! params.genome_accessions.isEmpty() ) { 
    log.info "Will download test data for SRA accessions: ${params.sra_accessions} and genome accessions: ${params.genome_accessions} using the TEST subworkflow for use in pipeline testing and development."
    TEST(
      ch_sra_accessions.ifEmpty([]), 
      ch_genome_accessions.ifEmpty([])
    )
    ch_reads    = ch_reads.mix(TEST.out.fastq)
    ch_fastas   = ch_fastas.mix(TEST.out.fasta)
    ch_versions = TEST.out.versions
  }

  log.info """
------------------------------------------------------

Initializing Complete

------------------------------------------------------------------------------------------------------------

"""

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

if ( ! params.skip_extras ) {
  workflow.onComplete {
    log.info """------------------------------------------------------

INITIALIZE subworkflow completed at: $workflow.complete

------------------------------------------------------
"""
  }
}
