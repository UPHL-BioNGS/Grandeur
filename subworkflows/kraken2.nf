include { kraken2_fastq as fastqs }  from '../modules/kraken2' addParams(params)
include { kraken2_fasta as contigs } from '../modules/kraken2' addParams(params)

workflow kraken2 {
    take:
        ch_reads
        ch_contigs
        ch_kraken2_db
  
    main:
        fastqs(ch_reads.combine(ch_kraken2_db))
        contigs(ch_contigs.combine(ch_kraken2_db))

    fastqs.out.results
        .mix(contigs.out.results)
        .map { it -> it [1] }
        .collectFile(
            storeDir: "${params.outdir}/kraken2/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "kraken2_summary.csv")
        .set { summary }

    emit:
        for_species = fastqs.out.results.mix(contigs.out.results)
        for_summary = summary
        for_multiqc = fastqs.out.for_multiqc.mix(contigs.out.for_multiqc)
}
