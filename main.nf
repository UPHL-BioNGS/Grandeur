#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UPHL-BioNGS/Grandeur
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UPHL-BioNGS/Grandeur
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INITIALIZE } from './subworkflows/local/initialize'
include { GRANDEUR   } from './workflows/grandeur'

include { paramsHelp; validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    main:
    //
    // HELP: Run help message and exit
    //

    if (params.help) {
        // We pass the usage string as a named parameter so it doesn't look for it in the schema keys
        def help_string = paramsHelp(usage: "nextflow run UPHL-BioNGS/Grandeur -profile docker --sample_sheet samplesheet.csv --outdir grandeur")
        log.info help_string
        exit 0
    }

    //
    // SUBWORKFLOW: Initialize files and tasks
    //
    INITIALIZE ()

    //
    // WORKFLOW: Run main workflow
    //
    GRANDEUR (
        INITIALIZE.out.reads,
        INITIALIZE.out.fastas,
        INITIALIZE.out.reference_genomes,
        INITIALIZE.out.versions,
        INITIALIZE.out.genome_sizes,
        INITIALIZE.out.mash_db,
        INITIALIZE.out.kraken2_db,
        INITIALIZE.out.checkm2_db,
        INITIALIZE.out.sylph_db,
        INITIALIZE.out.dataset_script,
        INITIALIZE.out.evaluat_script,
        INITIALIZE.out.jsoncon_script,
        INITIALIZE.out.multiqc_script,
        INITIALIZE.out.summary_script,
        INITIALIZE.out.summfle_script,
        INITIALIZE.out.version_script
    )

    log.info """
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Grandeur WORKFLOW                                                  ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ All results will be saved to:                                      │
│  ${params.outdir.padRight(66)}│
│    ├── multiqc/multiqc_report.html                                 │
│    └── grandeur_summary.tsv                                        │
└────────────────────────────────────────────────────────────────────┘

Thank you for using Grandeur! The view really is great from up here.

Please remember to tell us about all issues at https://github.com/UPHL-BioNGS/Grandeur/issues
------------------------------------------------------------------------------------------------------------
"""
}
