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

    workflow.onComplete {
        log.info """
------------------------------------------------------------------------------------------------------------

GRANDEUR pipeline execution summary
-----------------------------------
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Status       : ${workflow.success ? 'SUCCESS' : 'FAILED'}
Exit status  : ${workflow.exitStatus ?: 'N/A'}
"""

    if (workflow.success) {
        log.info """
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Completed Successfully                                    ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ All results have been saved to:                                    │
│  ${params.outdir.padRight(66)}│"""
        if ( ! params.skip_extras ) {
            log.info """│    ├── multiqc/multiqc_report.html                                 │
│    └── grandeur_summary.tsv                                        │"""
        }
        
        log.info """└────────────────────────────────────────────────────────────────────┘

Thanks for using Grandeur! The view really is great from up here.
------------------------------------------------------------------------------------------------------------
"""
    } else {
        log.info """
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Failed                                                    ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ Error message:                                                     │
│ ${workflow.errorMessage ?: 'No specific error message provided.'}
│                                                                    │
│ Please check the .nextflow.log file for more detailed information. │
└────────────────────────────────────────────────────────────────────┘

------------------------------------------------------------------------------------------------------------
"""
    }
    }
}