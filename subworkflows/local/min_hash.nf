include { MASH_SKETCH as SKETCH }  from '../../modules/local/mash'
include { MASH_DIST   as DIST   }  from '../../modules/local/mash'

workflow MIN_HASH {
    take:
        ch_reads
        ch_fastas
        ch_mash_db

    main:
        ch_versions      = Channel.empty()

        SKETCH(ch_reads.mix(ch_fastas))

        SKETCH.out.summary
            .collectFile(
                storeDir: "${params.outdir}/mash/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "mash_err_summary.csv")
            .set { mash_err_summary }

        ch_versions = ch_versions.mix(SKETCH.out.versions.first())

        if (params.mash_db) {
            DIST(SKETCH.out.msh.filter({it[1].size() > 0 }).combine(ch_mash_db))
        } else {
            DIST(SKETCH.out.msh.filter({it[1].size() > 0 }).map{it -> tuple(it[0], it[1], null)})
        }

        DIST.out.results
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/mash/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "mash_summary.csv")
            .set { mash_summary }

        ch_versions = ch_versions.mix(DIST.out.versions.first())

    emit:
        for_summary = mash_summary.mix(mash_err_summary)
        for_flag    = DIST.out.results
        versions    = ch_versions
}
