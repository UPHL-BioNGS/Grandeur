include { CIRCULOCOV }     from '../../modules/local/circulocov'
include { FASTQC }         from '../../modules/local/fastqc'
include { MLST }           from '../../modules/local/mlst'
include { PLASMIDFINDER }  from '../../modules/local/plasmidfinder'
include { QUAST }          from '../../modules/local/quast'

workflow QUALITY_ASSESSMENT {
    take:
    ch_reads
    ch_contigs
    ch_reads_contigs
    summfle_script

    main:
    for_multiqc = Channel.empty()
    ch_versions = Channel.empty()
    ch_summary  = Channel.empty()
    ch_bams     = Channel.empty()

    // fastq files only
    if ( params.sample_sheet || params.reads || params.sra_accessions ) {
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

        CIRCULOCOV(ch_reads_contigs.filter{it[1]}.filter{it[2]})
        ch_versions = ch_versions.mix(CIRCULOCOV.out.versions.first())

        CIRCULOCOV.out.collect
            .collectFile(name: "circulocov_summary.tsv",
                keepHeader: true,
                sort: { file -> file.text },
                storeDir: "${params.outdir}/circulocov")
            .set{ circulocov_summary }

        ch_summary  = ch_summary.mix(circulocov_summary)
        ch_bams     = ch_bams.mix(CIRCULOCOV.out.contig_bam)
    }

    // contigs
    QUAST(ch_reads_contigs.filter{it[2]})
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

    emit:
    bams        = ch_bams
    for_summary = ch_summary.collect()
    for_multiqc = for_multiqc.mix(QUAST.out.for_multiqc).collect()
    versions    = ch_versions
}
