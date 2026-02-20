include { DATASETS_SUMMARY }  from '../../../modules/local/datasets_summary'
include { DATASETS_DOWNLOAD } from '../../../modules/local/datasets_download'
include { REFERENCES }        from '../../../modules/local/references'
include { SKANI }             from '../../../modules/local/skani'
include { SPECIES }           from '../../../modules/local/species'
include { SPESTIMATOR }       from '../../../modules/local/spestimator'


workflow AVERAGE_NUCLEOTIDE_IDENTITY {
    take:
        ch_contigs
        ch_reference_genomes
        ch_species
        dataset_script

    main:
        log.info "Running average nucleotide identity (ANI) analysis)."
        ch_versions = channel.empty()
        
        if ( params.current_datasets ) {
            log.info "Downloading reference genomes for species in the dataset from NCBI with DATASETS."
            log.info "This is a third-party API that is not controlled by the Grandeur developers, requires the workflow to have internet access, and may be slow or have issues."
            log.info "Reference genomes are identified using results from SPESTIMATOR and MASH, as well as SYLPH and KRAKEN2 (if their respective databases are provided)."
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

            ch_datasets_summary
                .subscribe { summaryFile ->
                    def genomeCount = summaryFile.countLines() - 1 
                    log.info "Successfully downloaded ${genomeCount} genomes from NCBI for ANI analysis."
                    if (genomeCount == 0) {
                        log.warn "No genomes were downloaded from NCBI for ANI analysis. This may be due to issues with the NCBI API, or because no reference genomes were identified for the species in the dataset. If you believe there should be reference genomes available, you can try running this workflow again, or you can provide your own reference genomes with 'params.reference_genomes'."
                    }
                }

        } else {
            log.info "Using local files and those provided with this workflow."
            log.info "FYI: To download additional references from NCBI, set 'params.current_datasets'."
            ch_datasets_summary = channel.empty()
        }

        REFERENCES()

        ch_reference_genomes
            .mix(REFERENCES.out.fastas.flatten())
            .unique()
            .collect()
            .map { it -> tuple([it])}
            .set{ch_genomes}

        // SKANI(ch_contigs.combine(ch_genomes))

        // SKANI.out.results
        //     .map { it -> it [1] }
        //     .collectFile(
        //         storeDir: "${params.outdir}/skani/",
        //         keepHeader: true,
        //         sort: { file -> file.text },
        //         name: "skani_summary.csv")
        //     .set { ch_skani_summary }

        // skani_len_summary = channel.empty()

        // ch_versions = ch_versions.mix(SKANI.out.versions.first())

    emit:
        //for_summary = ch_skani_summary.mix(ch_datasets_summary).mix(skani_len_summary)
        for_summary = channel.empty()
        top_hit     = channel.empty()
        versions    = ch_versions
}

workflow.onComplete {
  log.info "Inititalization completed at: $workflow.complete"
  log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}