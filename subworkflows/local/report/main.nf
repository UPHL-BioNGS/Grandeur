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

        log.info """

Creating final reports

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process           ┃ description                                                        ┃
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ VERSIONS          ┃ Custom process to convert versions.yml for MultiQC.                ┃
┃ MULTIQC           ┃ Creation of html summary file.                                     ┃
┃ SUMMARY           ┃ Custom process that summarizes all results in text format.         ┃ 
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

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
        //summary  = SUMMARY.out.extended_tsv
        summary  = channel.empty()
        versions = ch_versions
}

if ( ! params.skip_extras ) {
    workflow.onComplete {
        log.info """------------------------------------------------------

REPORT subworkflow completed at: $workflow.complete

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Subworkflow Output Files                              ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│   'params.outdir'                                     │
│    ├── multiqc                                        │
│    │   └── multiqc_report.html                        │
│    └── grandeur_summary.tsv                           │
└───────────────────────────────────────────────────────┘

------------------------------------------------------
"""
    }
}
