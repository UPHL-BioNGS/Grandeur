include { BLASTN                     } from '../../modules/local/blast'
include { BLOBTOOLS_CREATE as CREATE } from '../../modules/local/blobtools'
include { BLOBTOOLS_PLOT as PLOT     } from '../../modules/local/blobtools'   
include { BLOBTOOLS_VIEW as VIEW     } from '../../modules/local/blobtools'
                                                                    
workflow BLOBTOOLS {
  take:
    ch_contig_bams
    ch_blast_db
  
  main:
    ch_versions = Channel.empty()
    ch_contigs  = ch_contig_bams.filter{it[1]}.map{it -> tuple(it[0], it[1])}
  
    BLASTN(ch_contigs.combine(ch_blast_db))
    ch_versions = ch_versions.mix(BLASTN.out.versions.first())

    CREATE(ch_contig_bams.join(BLASTN.out.blastn, by: 0, failOnMismatch: false, remainder: false))
    ch_versions = ch_versions.mix(CREATE.out.versions.first())

    VIEW(CREATE.out.json)
    ch_versions = ch_versions.mix(VIEW.out.versions.first())

    PLOT(CREATE.out.json)
    ch_versions = ch_versions.mix(PLOT.out.versions.first())
  
    PLOT.out.collect
      .collectFile(
        storeDir: "${params.outdir}/blobtools/",
        keepHeader: true,
        sort: { file -> file.text },
        name: "blobtools_summary.txt")
      .set{ summary }

  emit:
    for_flag    = PLOT.out.results
    for_summary = summary
    versions    = ch_versions
}
