include { ENA_DOWNLOAD as DOWNLOAD_FASTQ       } from '../../../modules/local/ena_download'
include { DATASETS_DOWNLOAD as DOWNLOAD_GENOME } from '../../../modules/local/datasets_download'

workflow TEST {
    take:
    ch_sra_accessions
    ch_genome_accessions

    main:

    log.info """

Downloading files from external databases.

More information can be found at https://github.com/UPHL-BioNGS/Grandeur/wiki/test.

Relevant params and their values:
- 'params.sra_accessions' : ${params.sra_accessions}
    - List of SRA accessions to download from the ENA
- 'params.genome_accessions' : ${params.genome_accessions}
    - List of genome accessions to download from NCBI genomes

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process           ┃ description                                                        ┃
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ DOWNLOAD_FASTQ    ┃ Downloads FASTQ files from ENA using enaDataGet                    ┃
┃ DOWNLOAD_GENOME   ┃ Downloads FASTA files from NCBI using DATASETS                     ┃
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

This subworkflow is dependant on third-part API which is out of the control by the 
Grandeur developers. These processes require internet access, and may be slow or have 
issues.

"""

    ch_versions = channel.empty()

    if ( ! params.sra_accessions.isEmpty() ) {
        DOWNLOAD_FASTQ(ch_sra_accessions.filter({it[0]}))
        ch_versions = ch_versions.mix(DOWNLOAD_FASTQ.out.versions.first())

        DOWNLOAD_FASTQ.out.fastq
            .map { it ->
                def meta = [id:it[0]] 
                tuple( meta, [file(it[1][0]), file(it[1][1])])
            }
            .set { ch_fastq }
    } else {
        ch_fastq = channel.empty()
    }

    if ( ! params.genome_accessions.isEmpty() ) {
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
        ch_fasta    = channel.empty()
    }

    emit:
    fastq    = ch_fastq
    fasta    = ch_fasta
    versions = ch_versions
}

if ( ! params.sra_accessions.isEmpty()  || ! params.genome_accessions.isEmpty() ) { 
    workflow.onComplete {
        log.info """------------------------------------------------------

TEST subworkflow completed at: $workflow.complete

------------------------------------------------------
"""
    }
}
