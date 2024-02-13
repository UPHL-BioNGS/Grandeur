include { mash_sketch_fastq }  from '../modules/local/mash' addParams(params)
include { mash_sketch_fasta }  from '../modules/local/mash' addParams(params)
include { mash_dist }          from '../modules/local/mash' addParams(params)
include { mash_err }           from '../modules/local/local' addParams(params)
//include { mash_screen } from '../modules/local/mash' addParams(params)

workflow min_hash {
    take:
        ch_reads
        ch_fastas
        ch_mash_db
  
    main:
        ch_mash_sketches = Channel.empty()
        ch_versions      = Channel.empty()

        if ( params.sample_sheet || params.reads ) {
            mash_sketch_fastq(ch_reads)

            // TODO : test mash_err 
            mash_err(mash_sketch_fastq.out.err)

            ch_mash_sketches = ch_mash_sketches.mix(mash_sketch_fastq.out.msh.filter({it[1].size() > 0 }))

            mash_err.out.summary
                .collectFile(
                    storeDir: "${params.outdir}/mash/",
                    keepHeader: true,
                    sort: { file -> file.text },
                    name: "mash_err_summary.csv")
                .set { mash_err_summary }

            ch_versions = ch_versions.mix(mash_sketch_fastq.out.versions.first())
        } else {
            mash_err_summary = Channel.empty()
        }

        if ( params.fastas || params.fasta_list ) {
            mash_sketch_fasta(ch_fastas)
            ch_mash_sketches = ch_mash_sketches.mix(mash_sketch_fasta.out.msh.filter({it[1].size() > 0 }))
            ch_versions      = ch_versions.mix(mash_sketch_fasta.out.versions.first())
        }

        if (params.mash_db) {
            mash_dist(ch_mash_sketches.combine(ch_mash_db))
        } else {
            mash_dist(ch_mash_sketches.map{it -> tuple(it[0], it[1], null)})
        }

        mash_dist.out.results
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/mash/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "mash_summary.csv")
            .set { mash_summary }

        ch_versions = ch_versions.mix(mash_dist.out.versions.first())

    emit:
        for_summary = mash_summary.mix(mash_err_summary)
        for_flag    = mash_dist.out.results
        versions    = ch_versions
}
