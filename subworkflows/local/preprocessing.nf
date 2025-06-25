/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP }   from '../../modules/local/fastp'

workflow PREPROCESSING {

    take:
    ch_reads

    main:
    ch_versions = Channel.empty()
    FASTP(ch_reads)
    
    FASTP.out.fastp_results
      .filter ({ it[2] as int >= params.minimum_reads })
      .map { it -> 
        tuple (it[0], it[1])
      }
      .set{ reads_check }

    ch_cleaned_reads = FASTP.out.fastp_results
    ch_versions      = ch_versions.mix(FASTP.out.versions.first())
    ch_multiqc       = FASTP.out.fastp_files

    emit:
    // TO-DO: FASTP.out.fastp is a subset of FASTP.out.fastp_results
    reads_check
    ch_cleaned_reads
    versions          = ch_versions
    for_multiqc       = ch_multiqc

}
