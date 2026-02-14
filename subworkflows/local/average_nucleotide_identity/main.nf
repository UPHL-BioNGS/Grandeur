include { DATASETS_SUMMARY }  from '../../../modules/local/datasets_summary'
include { DATASETS_DOWNLOAD } from '../../../modules/local/datasets_download'
include { REFERENCES }        from '../../../modules/local/references'
include { SKANI }             from '../../../modules/local/skani'
include { SPECIES }           from '../../../modules/local/species'
include { SPESTIMATOR }       from '../../../modules/local/spestimator'
include { SYLPH}              from '../../../modules/local/sylph' 

workflow AVERAGE_NUCLEOTIDE_IDENTITY {
    take:
        ch_species
        ch_contigs
        ch_reference_genomes
        ch_sylph_db
        dataset_script

    main:
        ch_versions = channel.empty()
        if ( params.current_datasets ) {
            SPECIES(ch_species)

            SPECIES.out.species
                .splitText()
                .map{ it -> it.trim()}
                .set{ ch_species_list }

            DATASETS_SUMMARY(ch_species_list.combine(dataset_script))
            DATASETS_DOWNLOAD(DATASETS_SUMMARY.out.genomes.collect())

            ch_reference_genomes = ch_reference_genomes.mix(DATASETS_DOWNLOAD.out.genomes.flatten())

            ch_versions = ch_versions.mix(DATASETS_SUMMARY.out.versions.first()).mix(DATASETS_DOWNLOAD.out.versions)

            DATASETS_SUMMARY.out.genomes
                .collectFile(
                    storeDir: "${params.outdir}/datasets/",
                    keepHeader: true,
                    sort: { file -> file.text },
                    name: "datasets_summary.csv")
                .set { ch_datasets_summary }

        } else {
            ch_datasets_summary = channel.empty()
        }

        if ( params.sylph_db ) {
            SYLPH(ch_contigs.combine(ch_sylph_db))
        }

        REFERENCES()

        ch_reference_genomes
            .mix(REFERENCES.out.fastas.flatten())
            .unique()
            .collect()
            .map { it -> tuple([it])}
            .set{ch_genomes}

        SKANI(ch_contigs.combine(ch_genomes))

        SKANI.out.results
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/skani/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "skani_summary.csv")
            .set { summary }

        SKANI.out.top_len
            .collectFile(
                keepHeader: true,
                name: "skani_top_len.csv")
            .set { skani_len_summary }

        ch_versions = ch_versions.mix(SKANI.out.versions.first())

    emit:
        for_flag    = SKANI.out.results
        for_summary = summary.mix(ch_datasets_summary).mix(skani_len_summary)
        top_hit     = channel.empty()
        versions    = ch_versions
}
