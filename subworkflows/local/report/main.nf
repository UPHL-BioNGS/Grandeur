include { MQC_PREP } from '../../../modules/local/mqc_prep'
include { MULTIQC }  from '../../../modules/local/multiqc'
include { SUMMARY }  from '../../../modules/local/summary'
include { VERSIONS } from '../../../modules/local/versions'

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

        ch_reads
            .mix(ch_fastas)
            .map { meta, files -> 
               def sample = meta.id
               def file1 = files[0].name
               def file2 = files[1] ? files[1].name : null
               def version = "${workflow.manifest.version}"
               return "${sample},${file1},${file2},${version}"
            }
            .collectFile(
                name: "input_files.txt",
                newLine: true
                )
            .set { ch_names }

        SUMMARY(for_summary.mix(ch_names).mix(MULTIQC.out.data_folder).collect())

    emit:
        summary  = SUMMARY.out.extended_tsv
        versions = ch_versions
}

if (! params.skip_extras ) {
    workflow.onComplete {
        log.info "Report workflow completed at: $workflow.complete"
        log.info "MultiQC report can be found at ${params.outdir}/multiqc/multiqc_report.html"
        log.info "Summary can be found at ${params.outdir}/grandeur_summary.tsv"
        log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }
}