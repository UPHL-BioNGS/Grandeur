include { datasets_summary }    from '../modules/datasets'  addParams(params)
include { datasets_download }   from '../modules/datasets'  addParams(params)
include { fastani }             from '../modules/fastani'   addParams(params)
include { species }             from '../modules/grandeur'  addParams(params)
include { decompression }       from '../modules/grandeur'  addParams(params)

workflow average_nucleotide_identity {
    take:
        ch_species
        ch_contigs
        ch_static_fastani_genomes
  
    main:

        if ( params.current_datasets ) {
            species(ch_species)

            species.out.species
                .splitText()
                .map(it -> it.trim())
                .set{ ch_species_list }

            datasets_summary(ch_species_list)
            datasets_download(datasets_summary.out.genomes.collect())

            ch_fastani_db = datasets_download.out.tar
        } else {
            ch_fastani_db = ch_static_fastani_genomes
        }

        decompression(ch_fastani_db)
        ch_contigs.combine(decompression.out.decompressed).view()
        fastani(ch_contigs.combine(decompression.out.decompressed))

   fastani.out.results
       .collectFile(
           storeDir: "${params.outdir}/fastani/",
           keepHeader: true,
           sort: { file -> file.text },
           name: "fastani_summary.csv")
       .set { summary }

   emit:
       for_summary = summary
}
