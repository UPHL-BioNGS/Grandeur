//
// Subworkflow with functionality specific to the nf-core/grandeur pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TEST } from "../../../subworkflows/local/test"

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_fastas    = Channel.empty()
    ch_versions  = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    // Create Channels for Scripts
    // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

    dataset_script = Channel.fromPath(workflow.projectDir + "/bin/datasets_download.py", type: "file")
    evaluat_script = Channel.fromPath(workflow.projectDir + "/bin/evaluate.py",          type: "file")
    jsoncon_script = Channel.fromPath(workflow.projectDir + "/bin/json_convert.py",      type: "file")
    multiqc_script = Channel.fromPath(workflow.projectDir + "/bin/for_multiqc.py",       type: "file")
    summary_script = Channel.fromPath(workflow.projectDir + "/bin/summary.py",           type: "file")
    summfle_script = Channel.fromPath(workflow.projectDir + "/bin/summary_file.py",      type: "file")
    version_script = Channel.fromPath(workflow.projectDir + "/bin/versions.py",          type: "file")


    // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    // Create Channels for User Inputs
    // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

    if (params.input){
        // using a sample sheet with the column header of 'sample,fastq_1,fastq_2'
        Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }
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

    if (params.fastas) {
        // getting fasta from a file
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
    }

    // Getting accession for downloading

    // from SRA
    ch_sra_accessions   = Channel.from( params.sra_accessions )

    // from genomes
    ch_genome_accessions = Channel.from( params.genome_accessions)

    // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    // Create Channels for Database Files
    // TODO: use error() instead of println()
    // ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

    // Getting the file with genome sizes of common organisms for fastqcscan. The End User can use their own file and set with a param
    Channel.fromPath(params.genome_sizes, type: "file").ifEmpty{
      println("The genome sizes file for this workflow are missing!! ${params.genome_sizes}")
      exit 1}
    .set { ch_genome_sizes }

    // Getting the database for blast/blobtools
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

    // User Supplied Fastani Reference Genomes
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

    // TODO: set ch_reads to reads obtained from ch_samplesheet, ch_sra_accessions, and ch_genome_accessions
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

