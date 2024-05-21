include { blastn }            from '../modules/local/blast'      addParams(params)
include { blobtools_create }  from '../modules/local/blobtools'  addParams(params)
include { blobtools_plot }    from '../modules/local/blobtools'  addParams(params)   
include { blobtools_view }    from '../modules/local/blobtools'  addParams(params)
                                                                    
workflow blobtools {
  take:
    ch_contig_bams
    ch_blast_db
  
  main:
    ch_contigs = ch_contig_bams.filter{it[1]}.map{it -> tuple(it[0], it[1])}
  
    blastn(ch_contigs.combine(ch_blast_db))
    blobtools_create(ch_contig_bams.join(blastn.out.blastn, by: 0, failOnMismatch: false, remainder: false))
    blobtools_view(blobtools_create.out.json)
    blobtools_plot(blobtools_create.out.json)
  
    blobtools_plot.out.collect
      .collectFile(
        storeDir: "${params.outdir}/blobtools/",
        keepHeader: true,
        sort: { file -> file.text },
        name: "blobtools_summary.txt")
      .set{ summary }

  emit:
    for_flag    = blobtools_plot.out.results
    for_summary = summary
    versions    = blastn.out.versions.first().mix(blobtools_create.out.versions.first()).mix(blobtools_view.out.versions.first()).mix(blobtools_plot.out.versions.first())
}
