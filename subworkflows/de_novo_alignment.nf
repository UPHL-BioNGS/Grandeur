include { fastp }   from '../modules/local/fastp'  addParams(params)
include { bbduk }   from '../modules/local/bbduk'  addParams(params)
include { spades }  from '../modules/local/spades' addParams(params)

workflow de_novo_alignment {
  take: 
    reads
  
  main:
    bbduk(reads)
    fastp(bbduk.out.fastq)

    fastp.out.fastq.view()

    fastp.out.fastq
      .map{it -> tuple(it , it[1][0].splitFastq( limit: params.minimum_reads , file: true) | count)}
      .view{ "From fastp : $it" }

    spades(fastp.out.fastq.filter{it[1][0].countFastq() >= params.minimum_reads})

  emit:
    // for downstream analyses
    clean_reads = fastp.out.fastq
    contigs     = spades.out.contigs

    // for multiqc
    for_multiqc = fastp.out.fastp_files.mix(bbduk.out.stats)
    versions    = bbduk.out.versions.mix(fastp.out.versions).mix(spades.out.versions)
}
