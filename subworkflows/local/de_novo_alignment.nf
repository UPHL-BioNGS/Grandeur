include { FASTP }   from '../../modules/local/fastp'
include { SPADES }  from '../../modules/local/spades'

workflow DE_NOVO_ALIGNMENT {
  take: 
    reads
  
  main:
    ch_versions = Channel.empty()

    FASTP(reads)

    FASTP.out.fastp_results
      .filter ({ it[2] as int >= params.minimum_reads })
      .map { it -> 
        tuple (it[0], it[1])
      }
      .set{ read_check }

    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    SPADES(read_check)

    ch_versions = ch_versions.mix(SPADES.out.versions.first())

  emit:
    // for downstream analyses
    reads_contigs = SPADES.out.reads_contigs
    clean_reads   = FASTP.out.fastq
    contigs       = SPADES.out.contigs.filter{it[1] != null}

    // for multiqc
    for_multiqc = FASTP.out.fastp_files
    versions    = ch_versions
}
