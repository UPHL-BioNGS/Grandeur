include { FASTP }   from '../../modules/local/fastp'
include { SPADES }  from '../../modules/local/spades'

workflow DE_NOVO_ALIGNMENT {
  take: 
    reads
  
  main:
    FASTP(reads)

    FASTP.out.fastp_results
      .filter ({ it[2] as int >= params.minimum_reads })
      .map { it -> 
        tuple (it[0], it[1])
      }
      .set{ read_check }

    SPADES(read_check)

  emit:
    // for downstream analyses
    reads_contigs = SPADES.out.reads_contigs
    clean_reads   = FASTP.out.fastq
    contigs       = SPADES.out.contigs.filter{it[1] != null}

    // for multiqc
    for_multiqc = FASTP.out.fastp_files
    versions    = FASTP.out.versions.first().mix(SPADES.out.versions.first())
}
