include { bbmap }             from '../modules/bbmap'      addParams(outdir: params.outdir)
include { blastn }            from '../modules/blast'      addParams(blast_db_type: params.blast_db_type, blastn_options: params.blastn_options )
include { blobtools_create }  from '../modules/blobtools'  addParams(blobtools_create_options: params.blobtools_create_options)
include { blobtools_plot }    from '../modules/blobtools'  addParams(blobtools_plot_options: params.blobtools_plot_options)    
include { blobtools_view }    from '../modules/blobtools'  addParams(blobtools_view_options: params.blobtools_view_options)
                                                                    
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
  
    blobtools_plot.out.results
      .collectFile(
        storeDir: "${params.outdir}/blobtools/",
        keepHeader: true,
        sort: true,
        name: "blobtools_species.txt")

  emit:
    species = blobtools_plot.out.species
    perc    = blobtools_plot.out.perc
}
