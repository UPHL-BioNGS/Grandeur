include { MASH_SKETCH as SKETCH }  from '../../modules/local/mash'
include { MASH_DIST   as DIST   }  from '../../modules/local/mash'
include { MASH_ERR    as ERR    }  from '../../modules/local/local'
//include { mash_screen } from '../modules/local/mash' addParams(params)

workflow MIN_HASH {
    take:
        ch_reads
        ch_fastas
        ch_mash_db

    main:
        ch_mash_sketches = Channel.empty()
        ch_versions      = Channel.empty()

        if ( params.sample_sheet || params.reads || params.sra_accessions ) {
            SKETCH(ch_reads)

            ERR(SKETCH.out.err)

            ch_mash_sketches = ch_mash_sketches.mix(SKETCH.out.msh.filter({it[1].size() > 0 }))

            ERR.out.summary
                .collectFile(
                    storeDir: "${params.outdir}/mash/",
                    keepHeader: true,
                    sort: { file -> file.text },
                    name: "mash_err_summary.csv")
                .set { mash_err_summary }

            ch_versions = ch_versions.mix(SKETCH.out.versions.first())
        } else {
            mash_err_summary = Channel.empty()
        }

        if ( params.fastas || params.fasta_list ) {
            SKETCH(ch_fastas)
            ch_mash_sketches = ch_mash_sketches.mix(SKETCH.out.msh.filter({it[1].size() > 0 }))
            ch_versions      = ch_versions.mix(SKETCH.out.versions.first())
        }

        if (params.mash_db) {
            DIST(ch_mash_sketches.combine(ch_mash_db))
        } else {
            DIST(ch_mash_sketches.map{it -> tuple(it[0], it[1], null)})
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
