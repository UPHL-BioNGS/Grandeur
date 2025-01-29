include { KRAKEN2 }  from '../../modules/local/kraken2'

workflow KMER_TAXONOMIC_CLASSIFICATION {
    take:
        ch_reads
        ch_kraken2_db

    main:
    KRAKEN2(ch_reads.combine(ch_kraken2_db))

    KRAKEN2.out.results
        .map { it -> it [1] }
        .collectFile(
            storeDir: "${params.outdir}/kraken2/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "kraken2_summary.csv")
        .set { summary }

    ch_versions = KRAKEN2.out.versions.first()

    emit:
        for_flag    = KRAKEN2.out.results
        for_summary = summary
        for_multiqc = KRAKEN2.out.for_multiqc
        versions    = ch_versions
}
