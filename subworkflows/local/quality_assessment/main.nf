include { AMRFINDER }      from '../../../modules/local/amrfinder'
include { CHECKM2 }         from '../../../modules/local/checkm2'
include { FASTQC }         from '../../../modules/local/fastqc'
include { MLST }           from '../../../modules/local/mlst'
include { PLASMIDFINDER }  from '../../../modules/local/plasmidfinder'
include { QUAST }          from '../../../modules/local/quast'

workflow QUALITY_ASSESSMENT {
    take:
    ch_reads
    ch_contigs
    ch_reads_contigs
    ch_checkm2_db
    summfle_script

    main:
    for_multiqc = Channel.empty()
    ch_versions = Channel.empty()
    ch_summary  = Channel.empty()
    ch_bams     = Channel.empty()

    log.info "Running quality assessment on the reads and assemblies. This workflow will perform quality control on the reads with FastQC, but the remaining processes of QUAST, CHECKM2, AMRFINDER, and PLASMIDFINDER will be run on generated assemblies as well as those specified with an input file designated with 'params.fasta_list'."

    // fastq files only, so hidden if only fasta files are provided
    if ( params.sample_sheet || params.reads || params.sra_accessions ) {
        log.info "Running quality assessment on the reads with FastQC. This will be performed on all read files provided, including those specified in the sample sheet, those provided with 'params.reads', and those downloaded from SRA with 'params.sra_accessions'."
        FASTQC(ch_reads)
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        for_multiqc = for_multiqc.mix(FASTQC.out.for_multiqc)

        FASTQC.out.collect
            .collectFile(name: "fastqc_summary.csv",
                keepHeader: true,
                sort: { file -> file.text },
                storeDir: "${params.outdir}/fastqc")
            .set{ fastqc_summary }

        ch_summary = ch_summary.mix(fastqc_summary)

    }

    AMRFINDER(ch_organism)

    AMRFINDER.out.collect
      .collectFile(name: 'amrfinderplus.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/amrfinder")
      .set{ amrfinderplus_summary }

    ch_summary  = ch_summary.mix(amrfinderplus_summary)
    ch_versions = ch_versions.mix(AMRFINDER.out.versions.first())

    QUAST(ch_reads_contigs)
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    QUAST.out.collect
        .collectFile(name: "quast_report.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/quast")
        .set{ quast_summary }

    ch_summary = ch_summary.mix(quast_summary)

    QUAST.out.collect_contig
        .collectFile(name: "quast_contig_report.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/quast")
        .set{ quast_contig_summary }
    ch_summary = ch_summary.mix(quast_contig_summary)

    MLST(ch_contigs.combine(summfle_script))
    ch_versions = ch_versions.mix(MLST.out.versions.first())

    MLST.out.collect
        .collectFile(name: "mlst_summary.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/mlst")
        .set{ mlst_summary }
    ch_summary = ch_summary.mix(mlst_summary)

    PLASMIDFINDER(ch_contigs.combine(summfle_script))
    ch_versions = ch_versions.mix(PLASMIDFINDER.out.versions.first())

    PLASMIDFINDER.out.collect
        .collectFile(name: "plasmidfinder_result.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/plasmidfinder")
        .set{ plasmidfinder_summary }
    ch_summary = ch_summary.mix(plasmidfinder_summary)

    if (params.checkm2_db) {    
        CHECKM2(ch_contigs.combine(ch_checkm2_db))
        ch_versions = ch_versions.mix(CHECKM2.out.versions.first())

        CHECKM2.out.collect
            .collectFile(name: "checkm2_summary.tsv",
                keepHeader: true,
                sort: { file -> file.text },
                storeDir: "${params.outdir}/checkm2")
            .set{ checkm2_summary }
        ch_summary = ch_summary.mix(checkm2_summary)
    } 

    emit:
    bams        = ch_bams
    for_summary = ch_summary.collect()
    for_multiqc = for_multiqc.mix(QUAST.out.for_multiqc).collect()
    versions    = ch_versions
}

workflow.onComplete {
  log.info "Inititalization completed at: $workflow.complete"
  log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}