include { CONCAT_REPORTS } from '../../../modules/local/concat_reports'
include { KRAKEN2     } from '../../../modules/local/kraken2'
include { MASH_DIST   } from '../../../modules/local/mashdist'
include { MASH_SCREEN } from '../../../modules/local/mashscreen'
include { SYLPH       } from '../../../modules/local/sylph'

workflow TAXONOMIC_PROFILING {
    take:
    ch_reads
    ch_fastas
    _ch_assemblies
    ch_kraken2_db
    ch_mash_db
    ch_sylph_db

    main:

        log.info """

This subworkflow is designed to take raw reads and/or assembled FASTA files and quickly 
identify the organisms present in the samples. It relies heavily on fast, 
k-mer/MinHash-based algorithms rather than computationally heavy alignments.

More information can be found at https://github.com/UPHL-BioNGS/Grandeur/wiki/taxonomic_profiling.

Relevant params and their values:
- 'params.kraken2_db' : ${params.kraken2_db}
    - Set to kraken2 directory
- 'params.mash_db' : ${params.mash_db}
    - Set to mash reference file
    - Used for both MASH_DIST and MASH_SCREEN
- 'params.sylph_db' : ${params.sylph_db}
    - Set to mash database file


┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process           ┃ description                                                        ┃
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ KRAKEN2           ┃ Uses KMERS to classify reads to taxa.                              ┃
┃ MASH_DIST         ┃ Uses MinHash sketches to rapidly estimate the distance between     ┃
┃                   ┃ genomic sequences.                                                 ┃
┃ MASH_SCREEN       ┃ Estimates containment of input files.                              ┃
┃ SYLPH             ┃ Uses a machine learning approach to classify reads to taxa.        ┃
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

    ch_versions = channel.empty()
    ch_summary  = channel.empty()
    ch_multiqc  = channel.empty()
    ch_species  = channel.empty()
    ch_concat   = channel.empty()

    if ( params.kraken2_db && ( params.sample_sheet || params.reads || params.sra_accessions )) {
        KRAKEN2(ch_reads.combine(ch_kraken2_db))
        ch_concat   = ch_concat.mix(KRAKEN2.out.results.map{ it -> it[1] }.collect().map {it -> [it, "kraken2_summary.csv","kraken2",true]})
        ch_versions = KRAKEN2.out.versions.first()
        ch_multiqc  = ch_multiqc.mix(KRAKEN2.out.for_multiqc)
    }

    if (params.mash_db) {
        MASH_DIST(ch_reads.mix(ch_fastas).filter{ it -> it }.combine(ch_mash_db))
        MASH_SCREEN(ch_reads.mix(ch_fastas).filter{ it -> it }.combine(ch_mash_db))
    } else {
        MASH_DIST(ch_reads.mix(ch_fastas).filter{ it -> it }.map{it -> tuple(it[0], it[1], null)})
        MASH_SCREEN(ch_reads.mix(ch_fastas).filter{ it -> it }.map{it -> tuple(it[0], it[1], null)})
    }
    ch_concat   = ch_concat.mix(MASH_DIST.out.results.collect().map {it -> [it, "mashdist_summary.csv","mash",false]})
    ch_concat   = ch_concat.mix(MASH_DIST.out.mash_err.collect().map {it -> [it, "mash_err_summary.csv","mash",false]})
    ch_concat   = ch_concat.mix(MASH_SCREEN.out.screen_results.collect().map {it -> [it, "skani_summary.csv","mash",false]})
    ch_versions = ch_versions.mix(MASH_DIST.out.versions.first())
    ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

    if (params.sylph_db) {
        SYLPH(ch_reads.mix(ch_fastas).filter{ it -> it }.combine(ch_sylph_db))
        ch_concat   = ch_concat.mix(SYLPH.out.results.collect().map {it -> [it, "sylph_summary.tsv","sylph",true]})
        ch_concat   = ch_concat.mix(SYLPH.out.for_download.collect().map {it -> [it, "sylph_download_summary.tsv","sylph",true]})
        ch_versions = ch_versions.mix(SYLPH.out.versions.first())
    }

    CONCAT_REPORTS(ch_concat)
    ch_species = ch_species.mix(CONCAT_REPORTS.out.summary.filter{ it -> it.contains("kraken") })
    ch_species = ch_species.mix(CONCAT_REPORTS.out.summary.filter{ it -> it.contains("mashscreen") })
    ch_species = ch_species.mix(CONCAT_REPORTS.out.summary.filter{ it -> it.contains("mashdist") })
    ch_species = ch_species.mix(CONCAT_REPORTS.out.summary.filter{ it -> it.contains("sylph_download_summary") })
    ch_summary = ch_summary.mix(CONCAT_REPORTS.out.summary)

    emit:
        for_ref_download = ch_species
        for_summary      = ch_summary
        for_multiqc      = ch_multiqc
        versions         = ch_versions
}

if ( ! params.skip_extras ) {
    workflow.onComplete {
        log.info """------------------------------------------------------

TAXONOMIC PROFILING subworkflow completed at: $workflow.complete

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Subworkflow Output Files                              ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│   ${params.outdir.padRight(52)}│"""
    if ( params.kraken2_db && ( params.sample_sheet || params.reads || params.sra_accessions )) {
        log.info """│    ├── kraken2                                        │
│    │   └── kraken2_summary.csv                        │"""
    }
    if (params.sylph_db ) {
        log.info """│    ├── sylph                                          │
│    │   └── sylph_summary.tsv                          │"""
        }

        log.info """│    └── mash                                           │
│        ├── mashdist_summary.csv                       │
│        └── mashscreen_summary.csv                     │
└───────────────────────────────────────────────────────┘

------------------------------------------------------
"""
    }
}
