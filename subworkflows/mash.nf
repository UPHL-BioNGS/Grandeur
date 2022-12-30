include { mash }  from '../modules/mash' addParams(params)

workflow mash {
    take:
        ch_reads
        ch_contigs
        ch_mash_db
  
    main:
        mash(ch_reads.mix(ch_contigs))

    mash.out.results
        .collectFile(
            storeDir: "${params.outdir}/mash/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "mash_summary.csv")
        .set { summary }

    emit:
        for_summary = summary
}
