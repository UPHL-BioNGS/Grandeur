include { FASTP }   from '../../../modules/local/fastp'
include { SPADES }  from '../../../modules/local/spades'

workflow DE_NOVO_ALIGNMENT {
  take: 
    reads
  
  main:
    ch_versions = Channel.empty()

    log.info "Running de novo assembly. This workflow will perform quality control on the reads with fastp, and then assemble the reads into contigs with SPAdes."
    log.info "The current threshold for the number of passed reads is ${params.minimum_reads}. Any samples with fewer than this will not be included in other steps."
    log.info "The minimum number of reads can be adjusted with 'params.minimum_reads'."

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

if ( params.sample_sheet || params.reads || params.sra_accessions ) {
  workflow.onComplete {
    log.info "Assembly completed at: $workflow.complete"
    log.info "Generated assemblies are at '${params.outdir}/contigs/'."
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
  }
}