include { FASTP }   from '../../modules/local/fastp'
include { SPADES }  from '../../modules/local/spades'

workflow DE_NOVO_ALIGNMENT {
  take: 
    reads
  
  main:
    ch_versions = Channel.empty()

    FASTP(reads)

    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    SPADES(FASTP.out.fastq)

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
