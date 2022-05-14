include { seqyclean } from '../modules/seqyclean' addParams(fastq_processes: params.fastq_processes, seqyclean_contaminant_file: params.seqyclean_contaminant_file, seqyclean_options: params.seqyclean_options)
include { spades }    from '../modules/spades'    addParams(fastq_processes: params.fastq_processes, spades_options: params.spades_options)

workflow de_novo_alignment {
  take: reads
  main:
    seqyclean(reads)
    spades(seqyclean.out.clean_reads)

    seqyclean.out.collect
      .collectFile(name: "SummaryStatistics.tsv",
      keepHeader: true,
      sort: true,
      storeDir: "${params.outdir}/seqyclean")

  emit:
    // for downstream analyses
    clean_reads               = seqyclean.out.clean_reads
    contigs                   = spades.out.contigs

    // for summary
    seqyclean_perc_kept       = seqyclean.out.perc_kept
    seqyclean_perc_pairskept  = seqyclean.out.pairskept

    // for multiqc
    seqyclean_multiqc         = seqyclean.out.for_multiqc
}
