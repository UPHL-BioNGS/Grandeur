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

if ( ! params.skip_extras ) {
  workflow.onComplete {
    log.info """------------------------------------------------------

GRANDEUR workflow completed at: $workflow.complete

Execution status: ${ workflow.success ? 'OK' : 'failed' }

------------------------------------------------------
"""
  }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
