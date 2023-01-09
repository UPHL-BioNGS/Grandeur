include { datasets_summary }    from '../modules/datasets'  addParams(params)
include { datasets_download }   from '../modules/datasets'  addParams(params)
include { decompression }       from '../modules/grandeur'  addParams(params)
include { fastani }             from '../modules/fastani'   addParams(params)
include { species }             from '../modules/grandeur'  addParams(params)

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
            ch_fastani_db = datasets_download.out.genomes

            datasets_summary.out.genomes
                .collectFile(
                   storeDir: "${params.outdir}/datasets/",
                    keepHeader: true,
                    sort: { file -> file.text },
                    name: "datasets_summary.csv")
                .set { datasets_genomes } 

        } else {
            decompression(ch_static_fastani_genomes)
            ch_fastani_db = decompression.out.decompressed
            datasets_genomes = Channel.empty()
        }

        fastani(ch_contigs.combine(ch_fastani_db))

        fastani.out.results
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/fastani/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "fastani_summary.csv")
            .set { summary }

   emit:
        for_species = fastani.out.results
        for_summary = summary
        for_size    = fastani.out.results.combine(datasets_genomes)
}
