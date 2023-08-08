include { kraken2 }  from '../modules/kraken2' addParams(params)

workflow kraken2 {
    take:
        ch_reads
        ch_kraken2_db
  
    main:
        kraken2(ch_reads.combine(ch_kraken2_db))

    kraken2.out.results
        .mix(contigs.out.results)
        .map { it -> it [1] }
        .collectFile(
            storeDir: "${params.outdir}/kraken2/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "kraken2_summary.csv")
        .set { summary }

    emit:
        for_flag    = kraken2.out.results.mix(contigs.out.results)
        for_summary = summary
        for_multiqc = kraken2.out.for_multiqc.mix(contigs.out.for_multiqc)
}
