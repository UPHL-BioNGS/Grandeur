include { fastp }     from '../modules/fastp'     addParams(fastp_options: params.fastp_options)
include { bbduk }     from '../modules/bbmap'     addParams(bbduk_options: params.bbduk_options)
include { spades }    from '../modules/spades'    addParams(spades_options: params.spades_options)

workflow de_novo_alignment {
  take: reads
  main:
    bbduk(reads)
    fastp(bbduk.out.fastq)

    fastp.out.fastq
      .join(fastp.out.fastp_results)
      .filter ({ it[2] as int >= params.minimum_reads })
      .map ( it -> tuple (it[0], it[1]))
      .set{ read_check }

    spades(read_check)

    // add fastqc here (since phoenix does it)

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
