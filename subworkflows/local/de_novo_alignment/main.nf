include { FASTP }   from '../../../modules/local/fastp'
include { SPADES }  from '../../../modules/local/spades'

workflow DE_NOVO_ALIGNMENT {
  take: 
    reads
  
  main:
    log.info """

Running de novo assembly.

More information can be found at https://github.com/UPHL-BioNGS/Grandeur/wiki/de_novo_alignment.

Relevant params and their values:
- 'params.minimum_reads' : ${params.minimum_reads}
    - Any samples with fewer than this will not be included in other steps.

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process           ┃ description                                                        ┃
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ FASTP             ┃ FASTQ files filtering                                              ┃
┃ SPADES            ┃ De novo assembly                                                   ┃
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

    ch_versions = channel.empty()

    FASTP(reads)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    SPADES(FASTP.out.fastq)
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    if ( params.sample_sheet || params.reads || params.sra_accessions ) {
      log.info """------------------------------------------------------

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ DE_NOVO_ALIGNMENT Output Files                        ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│   ${params.outdir.padRight(52)}│
│    └── contigs                                        │
│        └── *_contigs.fa                               │
└───────────────────────────────────────────────────────┘

------------------------------------------------------
"""
  }

  emit:
    // for downstream analyses
    reads_contigs = SPADES.out.reads_contigs
    clean_reads   = FASTP.out.fastq
    contigs       = SPADES.out.contigs.filter{it -> it[1] != null}

    // for multiqc
    for_multiqc = FASTP.out.fastp_files
    versions    = ch_versions
}


