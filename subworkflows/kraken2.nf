include { kraken2_fastq as FASTQ }  from '../modules/kraken2' addParams(params)
include { kraken2_fasta as CONTIG } from '../modules/kraken2' addParams(params)

workflow kraken2 {
    take:
        ch_reads
        ch_contigs
        ch_kraken2_db
  
    main:
        FASTQ(ch_reads.combine(ch_kraken2_db))
        CONTIG(ch_contigs.combine(ch_kraken2_db))

    FASTQ.out.results
        .mix(CONTIG.out.results)
        .collectFile(
            storeDir: "${params.outdir}/kraken2/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "kraken2_summary.csv")
        .set { summary }

    emit:
        for_summary = summary
        for_multiqc = FASTQ.out.for_multiqc.mix(CONTIG.out.for_multiqc)
}
