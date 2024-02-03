include { kraken2 }  from '../modules/local/kraken2' addParams(params)

workflow kmer_taxonomic_classification {
    take:
        ch_reads
        ch_kraken2_db
  
    main:
    kraken2(ch_reads.combine(ch_kraken2_db))

    kraken2.out.results
        .map { it -> it [1] }
        .collectFile(
            storeDir: "${params.outdir}/kraken2/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "kraken2_summary.csv")
        .set { summary }

    // TODO : get classified reads

    emit:
        for_flag    = kraken2.out.results
        for_summary = summary
        for_multiqc = kraken2.out.for_multiqc
        versions    = kraken2.out.versions
}
