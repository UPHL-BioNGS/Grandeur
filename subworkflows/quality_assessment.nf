include { circulocov }     from '../modules/local/circulocov'     addParams(params)
include { fastqc }         from '../modules/local/fastqc'         addParams(params)
include { mlst }           from '../modules/local/mlst'           addParams(params)
include { plasmidfinder }  from '../modules/local/plasmidfinder'  addParams(params)
include { quast }          from '../modules/local/quast'          addParams(params)

workflow quality_assessment {
  take:
    ch_reads
    ch_contigs
    summfle_script

  main:
    for_multiqc = Channel.empty()
    ch_versions = Channel.empty()
    ch_summary  = Channel.empty()
    ch_bams     = Channel.empty()

    // fastq files only
    if ( params.sample_sheet || params.reads || params.sra_accessions ) {
        fastqc(ch_reads)

        ch_reads
            .join(ch_contigs, by: 0, remainder: true)
            .filter {it[1]}
            .filter {it[2]}
            .set { for_circulocov }

        circulocov(for_circulocov)

        for_multiqc = for_multiqc.mix(fastqc.out.for_multiqc)

        fastqc.out.collect
            .collectFile(name: "fastqc_summary.csv",
                keepHeader: true,
                sort: { file -> file.text },
                storeDir: "${params.outdir}/fastqc")
            .set{ fastqc_summary }

        circulocov.out.collect
            .collectFile(name: "circulocov_summary.tsv",
                keepHeader: true,
                sort: { file -> file.text },
                storeDir: "${params.outdir}/circulocov")
            .set{ circulocov_summary }

        ch_summary  = ch_summary.mix(circulocov_summary).mix(fastqc_summary)
        ch_versions = ch_versions.mix(fastqc.out.versions.first()).mix(circulocov.out.versions.first())
        ch_bams     = ch_bams.mix(circulocov.out.bam)

    }

    // contigs
    quast(ch_contigs.join(ch_reads, by: 0, remainder: true ))
    mlst(ch_contigs.combine(summfle_script))    
    plasmidfinder(ch_contigs.combine(summfle_script))

    mlst.out.collect
        .collectFile(name: "mlst_summary.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/mlst")
        .set{ mlst_summary }

    plasmidfinder.out.collect
        .collectFile(name: "plasmidfinder_result.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/plasmidfinder")
        .set{ plasmidfinder_summary }

    quast.out.collect
        .collectFile(name: "quast_report.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/quast")
        .set{ quast_summary }

    quast.out.collect_contig
        .collectFile(name: "quast_contig_report.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/quast")
        .set{ quast_contig_summary }

    ch_summary
        .mix(mlst_summary)
        .mix(plasmidfinder_summary)
        .mix(quast_summary)
        .mix(quast_contig_summary)
        .set { for_summary }

    ch_versions
      .mix(mlst.out.versions.first())
      .mix(plasmidfinder.out.versions.first())
      .mix(quast.out.versions.first())
      .set { for_versions }

  emit:
    bams        = ch_bams
    for_summary = for_summary.collect()
    for_multiqc = for_multiqc.mix(quast.out.for_multiqc).collect()
    versions    = for_versions
}
