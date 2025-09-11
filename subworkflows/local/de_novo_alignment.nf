/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPADES }  from '../../modules/local/spades'

workflow DE_NOVO_ALIGNMENT {
    take:
    reads_check
    ch_versions

    main:
    
    SPADES(reads_check)

    reads_contigs = SPADES.out.reads_contigs
    ch_contigs  = SPADES.out.contigs.filter{it[1] != null}
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    emit:

    reads_contigs = reads_contigs
    contigs       = ch_contigs
    versions      = ch_versions

}
