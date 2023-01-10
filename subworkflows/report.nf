include { names }   from '../modules/summary' addParams(params)
include { multiqc } from '../modules/multiqc' addParams(params)
include { summary } from '../modules/summary' addParams(params)

workflow report {
    take:
        ch_reads
        ch_fastas
        for_multiqc
        for_summary
        fastp_reads
        phix_reads
  
    main:
        multiqc(for_multiqc)

        names(ch_reads.mix(ch_fastas).join(fastp_reads, by: 0, remainder: true).join(phix_reads, by: 0 , remainder: true))

        names.out.collect
        .collectFile(
           storeDir: "${params.outdir}/summary/",
           keepHeader: true,
           sort: { file -> file.text },
           name: "input_files.csv")
        .set { ch_names }

        summary(for_summary.mix(ch_names).collect())
}