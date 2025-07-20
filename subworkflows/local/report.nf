include { NAMES }    from '../../modules/local/local'
include { MQC_PREP } from '../../modules/local/local'
include { MULTIQC }  from '../../modules/local/multiqc'
include { SUMMARY }  from '../../modules/local/local'
include { VERSIONS } from '../../modules/local/multiqc'

workflow REPORT {
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

        VERSIONS(ch_collated_versions, version_script)

        MQC_PREP(for_multiqc.mix(for_summary).collect(), multiqc_script)

        MULTIQC(for_multiqc.mix(for_summary).mix(MQC_PREP.out.for_multiqc).mix(VERSIONS.out.for_multiqc).collect())

        NAMES(ch_reads.mix(ch_fastas))

        NAMES.out.collect
            .collectFile(
                keepHeader: true,
                sort: { file -> file.text },
                name: "input_files.csv")
            .set { ch_names }

        SUMMARY(for_summary.mix(ch_names).mix(MULTIQC.out.data_folder).collect())

    emit:
        summary  = SUMMARY.out.extended_tsv
        versions = ch_versions
}