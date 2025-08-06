include { DATASETS_SUMMARY }  from '../../modules/local/datasets'
include { DATASETS_DOWNLOAD } from '../../modules/local/datasets'
include { FASTANI }           from '../../modules/local/fastani'
include { REFERENCES }        from '../../modules/local/local'
include { SPECIES }           from '../../modules/local/local'

workflow AVERAGE_NUCLEOTIDE_IDENTITY {
    take:
        ch_species
        ch_contigs
        ch_fastani_ref
        dataset_script

    main:
        ch_versions = Channel.empty()
        if ( params.current_datasets ) {
            SPECIES(ch_species)

            SPECIES.out.species
                .splitText()
                .map{ it -> it.trim()}
                .set{ ch_species_list }

            DATASETS_SUMMARY(ch_species_list.combine(dataset_script))
            DATASETS_DOWNLOAD(DATASETS_SUMMARY.out.genomes.collect())

            ch_fastani_ref = ch_fastani_ref.mix(DATASETS_DOWNLOAD.out.genomes.flatten())

            ch_versions = ch_versions.mix(DATASETS_SUMMARY.out.versions.first()).mix(DATASETS_DOWNLOAD.out.versions)

            DATASETS_SUMMARY.out.genomes
                .collectFile(
                    storeDir: "${params.outdir}/datasets/",
                    keepHeader: true,
                    sort: { file -> file.text },
                    name: "datasets_summary.csv")
                .set { ch_datasets_summary }

        } else {
            ch_datasets_summary = Channel.empty()
        }

        REFERENCES()

        ch_fastani_ref
            .mix(REFERENCES.out.fastas.flatten())
            .unique()
            .collect()
            .map { it -> tuple([it])}
            .set{ch_fastani_genomes}

        FASTANI(ch_contigs.combine(ch_fastani_genomes))

        FASTANI.out.results
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/fastani/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "fastani_summary.csv")
            .set { summary }

        FASTANI.out.top_len
            .collectFile(
                keepHeader: true,
                name: "fastani_top_len.csv")
            .set { fastani_len_summary }

        ch_versions = ch_versions.mix(FASTANI.out.versions.first())

    emit:
        for_flag    = FASTANI.out.results
        for_summary = summary.mix(ch_datasets_summary).mix(fastani_len_summary)
        top_hit     = FASTANI.out.top_hit.map{ it -> tuple(it[0], it[1].baseName, it[1])}
        versions    = ch_versions
}
