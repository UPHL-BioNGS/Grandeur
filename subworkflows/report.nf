include { names }   from '../modules/grandeur' addParams(params)
include { multiqc } from '../modules/multiqc'  addParams(params)
include { summary } from '../modules/grandeur' addParams(params)

workflow report {
    take:
        ch_reads
        ch_fastas
        for_multiqc
        for_summary
  
    main:
        multiqc(for_multiqc.mix(for_summary).collect())

        names(ch_reads.mix(ch_fastas))

        names.out.collect
            .collectFile(
                storeDir: "${params.outdir}/summary/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "input_files.csv")
            .set { ch_names }

        summary(for_summary.mix(ch_names).mix(multiqc.out.data_folder).collect())
}