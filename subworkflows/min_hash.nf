include { mash_sketch_fastq }  from '../modules/local/mash' addParams(params)
include { mash_sketch_fasta }  from '../modules/local/mash' addParams(params)
include { mash_dist }   from '../modules/local/mash' addParams(params)
//include { mash_screen } from '../modules/local/mash' addParams(params)

workflow min_hash {
    take:
        ch_reads
        ch_fastas
        ch_mash_db
  
    main:
        mash_sketch_fastq(ch_reads)
        mash_sketch_fasta(ch_fastas)

        mash_sketch_fastq.out.msh
            .mix(mash_sketch_fasta.out.msh)
            .filter({it[1].size() > 0 })
            .set { ch_mash_sketches }

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
            .set { summary }

    emit:
        for_summary = summary
        for_flag    = mash_dist.out.results
        mash_err    = mash_sketch_fasta.out.err.mix(mash_sketch_fastq.out.err)
}
