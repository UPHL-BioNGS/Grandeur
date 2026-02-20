include { ENA_DOWNLOAD as DOWNLOAD_FASTQ       } from '../../../modules/local/ena_download'
include { DATASETS_DOWNLOAD as DOWNLOAD_GENOME } from '../../../modules/local/datasets_download'

workflow TEST {
    take:
    ch_sra_accessions
    ch_genome_accessions

    main:
    ch_versions = Channel.empty()

    if ( ! params.sra_accessions.isEmpty() ) {
        log.info "Downloading FASTQ files using enaDataGet for the following accessions: ${params.sra_accessions}."
        log.info "This is a third-party API that is not controlled by the Grandeur developers, requires the workflow to have internet access, and may be slow or have issues."
        DOWNLOAD_FASTQ(ch_sra_accessions.filter({it[0]}))
        ch_versions = ch_versions.mix(DOWNLOAD_FASTQ.out.versions.first())

        DOWNLOAD_FASTQ.out.fastq
            .map { it ->
                def meta = [id:it[0]] 
                tuple( meta, [file(it[1][0]), file(it[1][1])])
            }
            .set { ch_fastq }
    } else {
        ch_fastq = Channel.empty()
    }

    if ( ! params.genome_accessions.isEmpty() ) {
        log.info "Downloading FASTA files from NCBI for the following accessions: ${params.genome_accessions}."
        log.info "This is a third-party API that is not controlled by the Grandeur developers, requires the workflow to have internet access, and may be slow or have issues."
        DOWNLOAD_GENOME(ch_genome_accessions.collectFile(name: 'ids.csv', newLine: true))
        ch_versions = ch_versions.mix(DOWNLOAD_GENOME.out.versions.first())

        DOWNLOAD_GENOME.out.genomes
            .flatten()
            .map { it ->
                def meta = [id:it.baseName]
                tuple( meta, it)
            }
            .set { ch_fasta }
    } else {
        ch_fasta    = Channel.empty()
    }

    emit:
    fastq    = ch_fastq
    fasta    = ch_fasta
    versions = ch_versions
}

workflow.onComplete {
  log.info "Inititalization completed at: $workflow.complete"
  log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}