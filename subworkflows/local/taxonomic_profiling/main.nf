include { KRAKEN2  } from '../../../modules/local/kraken2'
include { MASH     } from '../../../modules/local/mash'
include { SYLPH    } from '../../../modules/local/sylph' 
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
        MASH(ch_reads.mix(ch_fastas).filter { it }.combine(ch_mash_db))
    } else {
        MASH(ch_reads.mix(ch_fastas).filter { it }.map{it -> tuple(it[0], it[1], null)})
    }

    MASH.out.results
        .map { it -> it [1] }
        .collectFile(
            storeDir: "${params.outdir}/mash/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "mash_summary.csv")
        .set { ch_mash_summary }

    ch_versions = ch_versions.mix(MASH.out.versions.first())
    ch_species  = ch_species.mix(ch_mash_summary)

    if (params.sylph_db) {
        SYLPH(ch_reads.mix(ch_fastas).filter { it }.combine(ch_sylph_db))

        SYLPH.out.tsv
        .map { it -> it [1] }
        .collectFile(
            storeDir: "${params.outdir}/sylph/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "sylph_summary.tsv")
        .set { ch_sylph_summary }

        ch_versions = ch_versions.mix(SYLPH.out.versions.first())
        ch_summary  = ch_summary.mix(ch_sylph_summary)
        //ch_multiqc  = ch_multiqc.mix(SYLPH.out.for_multiqc)
        //ch_species  = ch_species.mix(SYLPH.out.species)
    }


    emit:
        for_ref_download = ch_species
        for_summary      = ch_summary
        for_multiqc      = ch_multiqc
        versions         = ch_versions
}
