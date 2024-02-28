include { names }    from '../modules/local/local'   addParams(params)
include { mqc_prep } from '../modules/local/local'   addParams(params)
include { multiqc }  from '../modules/local/multiqc' addParams(params)
include { summary }  from '../modules/local/local'   addParams(params)
include { versions } from '../modules/local/multiqc' addParams(params)

workflow report {
    take:
        ch_reads
        ch_fastas
        for_multiqc
        for_summary
        ch_versions
        multiqc_script
        version_script
  
    main:
        ch_versions
            .collectFile(
                keepHeader: false,
                name: "versions.yml")
            .set { ch_collated_versions }

        versions(ch_collated_versions, version_script)

        mqc_prep(for_multiqc.mix(for_summary).collect(), multiqc_script)

        multiqc(for_multiqc.mix(for_summary).mix(mqc_prep.out.for_multiqc).mix(versions.out.for_multiqc).collect())

        //multiqc(for_multiqc.mix(for_summary).mix(mqc_prep.out.for_multiqc).collect())

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