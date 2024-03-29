include { datasets_summary }  from '../modules/local/datasets' addParams(params)
include { datasets_download } from '../modules/local/datasets' addParams(params)
include { fastani }           from '../modules/local/fastani'  addParams(params)
include { references }        from '../modules/local/local'    addParams(params)
include { species }           from '../modules/local/local'    addParams(params)

workflow average_nucleotide_identity {
    take:
        ch_species
        ch_contigs
        ch_fastani_ref
        dataset_script
  
    main:
        ch_versions = Channel.empty()
        if ( params.current_datasets ) {
            species(ch_species)

            species.out.species
                .splitText()
                .map(it -> it.trim())
                .set{ ch_species_list }

            datasets_summary(ch_species_list.combine(dataset_script))
            datasets_download(datasets_summary.out.genomes.collect())

            ch_fastani_ref = ch_fastani_ref.mix(datasets_download.out.genomes.flatten())

            ch_versions = ch_versions.mix(datasets_summary.out.versions.first()).mix(datasets_download.out.versions)

            datasets_summary.out.genomes
                .collectFile(
                   storeDir: "${params.outdir}/datasets/",
                    keepHeader: true,
                    sort: { file -> file.text },
                    name: "datasets_summary.csv")
                .set { datasets_summary }

        } else {
            datasets_summary = Channel.empty()
        }

        references()

        ch_fastani_ref
            .mix(references.out.fastas.flatten())
            .unique()
            .collect()
            .map { it -> tuple([it])}
            .set{ch_fastani_genomes}

        fastani(ch_contigs.combine(ch_fastani_genomes))

        fastani.out.results
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/fastani/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "fastani_summary.csv")
            .set { summary }

        fastani.out.top_len
            .collectFile(
                keepHeader: true,
                name: "fastani_top_len.csv")
            .set { fastani_len_summary }

        ch_versions = ch_versions.mix(fastani.out.versions.first())

    emit:
        for_flag    = fastani.out.results
        for_summary = summary.mix(datasets_summary).mix(fastani_len_summary)
        top_hit     = fastani.out.top_hit
        versions    = ch_versions
}
