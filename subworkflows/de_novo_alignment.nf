include { fastp }   from '../modules/local/fastp'  addParams(params)
include { bbduk }   from '../modules/local/bbduk'  addParams(params)
include { spades }  from '../modules/local/spades' addParams(params)

workflow de_novo_alignment {
  take: 
    reads
  
  main:
    bbduk(reads)
    fastp(bbduk.out.fastq)

    fastp.out.fastp_results
      .filter ({ it[2] as int >= params.minimum_reads })
      .map ( it -> tuple (it[0], it[1]))
      .set{ read_check }

    spades(read_check)

  emit:
    // for downstream analyses
    reads_contigs = spades.out.reads_contigs
    clean_reads   = fastp.out.fastq
    contigs       = spades.out.contigs.filter{it[1] != null}

    // for multiqc
    for_multiqc = fastp.out.fastp_files.mix(bbduk.out.stats)
    versions    = bbduk.out.versions.first().mix(fastp.out.versions.first()).mix(spades.out.versions.first())
}
