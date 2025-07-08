/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_grandeur_pipeline'

// GRANDEUR PIPELINE SUBWORKFLOWS
include { PREPROCESSING          } from '../subworkflows/local/preprocessing'
include { DE_NOVO_ALIGNMENT      } from '../subworkflows/local/de_novo_alignment'
include { QUALITY_ASSESSMENT     } from '../subworkflows/local/quality_assessment'
include { MIN_HASH               } from "../subworkflows/local/min_hash"
include { BLOBTOOLS              } from "../subworkflows/local/blobtools"
include { KMER_TAXONOMIC_CLASSIFICATION } from "../subworkflows/local/kmer_taxonomic_classification"

// GRANDEUR PIPELINE MODULES

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GRANDEUR {

    take:
    // ch_samplesheet // channel: samplesheet read in from --input
    ch_raw_reads
    ch_fastas
    ch_fastani_genomes
    ch_versions
    ch_genome_sizes
    ch_mash_db
    ch_kraken2_db
    ch_blast_db
    dataset_script
    evaluat_script
    jsoncon_script
    multiqc_script
    summary_script
    summfle_script
    version_script

    main:

    ch_for_multiqc   = Channel.empty()
    ch_for_summary   = ch_genome_sizes
    ch_versions      = Channel.empty()
    // ch_multiqc_files = Channel.empty()
    ch_for_flag      = Channel.empty()
    ch_reads_contigs = ch_fastas.map{it -> tuple(it[0], it[1], null)}

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        GRANDEUR PIPELINE LOGIC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if ( params.sample_sheet || params.reads || params.sra_accessions ) {
        PREPROCESSING (
            ch_raw_reads
        )

        reads_check        = PREPROCESSING.out.reads_check
        ch_clean_reads     = PREPROCESSING.out.ch_cleaned_reads
        ch_versions        = ch_versions.mix(PREPROCESSING.out.versions)

        ch_for_multiqc     = ch_for_multiqc.mix(PREPROCESSING.out.for_multiqc)

        DE_NOVO_ALIGNMENT (
            reads_check,
            ch_versions
        )

        // TODO: ch_contigs is a subset of ch_reads_contigs
        ch_contigs         = ch_fastas.mix(DE_NOVO_ALIGNMENT.out.contigs)
        ch_reads_contigs   = ch_reads_contigs.mix(DE_NOVO_ALIGNMENT.out.reads_contigs)
        ch_versions        = ch_versions.mix(DE_NOVO_ALIGNMENT.out.versions)

    } else {
        ch_contigs         = ch_fastas
        ch_cleaned_reads   = Channel.empty()
    }

    // getting a summary of everything
    if ( ! params.skip_extras ) {
        QUALITY_ASSESSMENT(
            ch_raw_reads,
            ch_contigs,
            ch_reads_contigs,
            summfle_script
        )

        ch_for_multiqc = ch_for_multiqc.mix(QUALITY_ASSESSMENT.out.for_multiqc)
        ch_for_summary = ch_for_summary.mix(QUALITY_ASSESSMENT.out.for_summary)
        ch_versions    = ch_versions.mix(QUALITY_ASSESSMENT.out.versions)

        // optional subworkflow blobtools (useful for interspecies contamination)
        if ( params.blast_db && ( params.sample_sheet || params.reads || params.sra_accessions )) {
            BLOBTOOLS(QUALITY_ASSESSMENT.out.bams, ch_blast_db )

            ch_for_summary = ch_for_summary.mix(BLOBTOOLS.out.for_summary)
            ch_for_flag    = ch_for_flag.mix(BLOBTOOLS.out.for_flag)
            ch_versions = ch_versions.mix(BLOBTOOLS.out.versions)
        }

        // optional subworkflow kraken2 (useful for interspecies contamination)
        if ( params.kraken2_db && ( params.sample_sheet || params.reads || params.sra_accessions )) {
            KMER_TAXONOMIC_CLASSIFICATION(ch_clean_reads, ch_kraken2_db )

            ch_for_multiqc = ch_for_multiqc.mix(KMER_TAXONOMIC_CLASSIFICATION.out.for_multiqc)
            ch_for_summary = ch_for_summary.mix(KMER_TAXONOMIC_CLASSIFICATION.out.for_summary)
            ch_for_flag    = ch_for_flag.mix(KMER_TAXONOMIC_CLASSIFICATION.out.for_flag)
            ch_versions    = ch_versions.mix(KMER_TAXONOMIC_CLASSIFICATION.out.versions)
        }

        // subworkflow mash for species determination
        MIN_HASH(ch_clean_reads, ch_fastas, ch_mash_db)
        ch_versions = ch_versions.mix(MIN_HASH.out.versions)
        ch_for_summary = ch_for_summary.mix(MIN_HASH.out.for_summary)
    }

}

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // Collate and save software versions
    //
    /*softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'grandeur_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
