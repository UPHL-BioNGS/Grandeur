include { bbmap }             from '../modules/bbmap'      addParams(params)
include { blastn }            from '../modules/blast'      addParams(params)
include { blobtools_create }  from '../modules/blobtools'  addParams(params)
include { blobtools_plot }    from '../modules/blobtools'  addParams(params)   
include { blobtools_view }    from '../modules/blobtools'  addParams(params)
                                                                    
workflow blobtools {
  take:
    clean_reads
    contigs
    blast_db
  
  main:
    bbmap(clean_reads.join(contigs, by: 0))
    blastn(contigs.combine(blast_db))
    blobtools_create(contigs.join(blastn.out.blastn, by: 0).join(bbmap.out.bam, by: 0))
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
    for_multiqc = bbmap.out.stats
}
