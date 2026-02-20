include { KRAKEN2  } from '../../../modules/local/kraken2'
include { MASH     } from '../../../modules/local/mash'
include { SYLPH    } from '../../../modules/local/sylph'
// still struggling with downloading the database
//include { SOURMASH } from '../../../modules/local/sourmash'

workflow TAXONOMIC_PROFILING {
    take:
    ch_reads
    ch_fastas
    ch_assemblies
    ch_kraken2_db
    ch_mash_db
    ch_sylph_db

    main:
    ch_versions = channel.empty()
    ch_summary  = channel.empty()
    ch_multiqc  = channel.empty()
    ch_species  = channel.empty()

    if ( params.kraken2_db && ( params.sample_sheet || params.reads || params.sra_accessions )) {
        log.info "KRAKEN2 uses KMERS to classify reads to taxa. The database used for classification can be adjusted with 'params.kraken2_db'."
        KRAKEN2(ch_reads.combine(ch_kraken2_db))

        KRAKEN2.out.results
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/kraken2/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "kraken2_summary.csv")
            .set { ch_kraken2_summary }

        ch_versions = KRAKEN2.out.versions.first()
        ch_multiqc  = ch_multiqc.mix(KRAKEN2.out.for_multiqc)
        ch_species  = ch_species.mix(ch_kraken2_summary)
    }

    if (params.mash_db) {
        log.info "MASH uses MinHash sketches to rapidly estimate the distance between genomic sequences. The database used for comparison can be adjusted with 'params.mash_db'."
        MASH(ch_reads.mix(ch_fastas).filter { it }.combine(ch_mash_db))
    } else {
        MASH(ch_reads.mix(ch_fastas).filter { it }.map{it -> tuple(it[0], it[1], null)})
    }

    MASH.out.results
        .collectFile(
            storeDir: file("${params.outdir}/mash/"),
            keepHeader: true,
            sort: { file -> file.text },
            name: "mashdist_summary.csv")
        .set { ch_mashdist_summary }

    MASH.out.screen_results
        .collectFile(
            storeDir: file("${params.outdir}/mash/"),
            keepHeader: true,
            sort: { file -> file.text },
            name: "mashscreen_summary.csv")
        .set { ch_mashscreen_summary }

    ch_versions = ch_versions.mix(MASH.out.versions.first())
    ch_species  = ch_species.mix(ch_mashdist_summary).mix(ch_mashscreen_summary)

    if (params.sylph_db) {
        log.info "SYLPH uses a machine learning approach to classify reads to taxa. The database used for classification can be adjusted with 'params.sylph_db'."
        SYLPH(ch_reads.mix(ch_fastas).filter { it }.combine(ch_sylph_db))

        SYLPH.out.tsv
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/sylph/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "sylph_summary.tsv")
            .set { ch_sylph_summary }

        SYLPH.out.for_download
            .collectFile(
                storeDir: "${params.outdir}/summary/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "sylph_download_summary.tsv")
            .set { ch_sylph_download_summary }

        ch_versions = ch_versions.mix(SYLPH.out.versions.first())
        ch_summary  = ch_summary.mix(ch_sylph_summary)
        ch_species  = ch_species.mix(ch_sylph_download_summary)
    }


    emit:
        for_ref_download = ch_species
        for_summary      = ch_summary
        for_multiqc      = ch_multiqc
        versions         = ch_versions
}

workflow.onComplete {
  log.info "Inititalization completed at: $workflow.complete"
  log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}