include { fastp }     from '../modules/fastp'     addParams(fastp_options: params.fastp_options)
include { bbduk }     from '../modules/bbmap'     addParams(bbduk_options: params.bbduk_options)
include { spades }    from '../modules/spades'    addParams(fastq_processes: params.fastq_processes, spades_options: params.spades_options)

workflow de_novo_alignment {
  take: reads
  main:
    bbduk(reads)
    fastp(bbduk.out.fastq)
    spades(fastp.out.fastq)

  emit:
    // for summary
    fastp_results             = fastp.out.fastp_results
    phix_reads                = bbduk.out.phix_reads

    // for downstream analyses
    clean_reads               = fastp.out.fastq
    contigs                   = spades.out.contigs

    // for multiqc
    fastp_multiqc             = fastp.out.fastp_files
    bbduk_multiqc             = bbduk.out.stats
}
