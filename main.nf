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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    main:
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


}

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
│  📁 ${params.outdir}                                               │
"""
        // Only point out the summary and MultiQC if they were actually generated
        if ( ! params.skip_extras ) {
            log.info """│                                                                    │
│ Key consolidated reports to check:                                 │
│  📄 ${params.outdir}/grandeur_summary.tsv                          │
│  📊 ${params.outdir}/multiqc/multiqc_report.html                   │"""
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
