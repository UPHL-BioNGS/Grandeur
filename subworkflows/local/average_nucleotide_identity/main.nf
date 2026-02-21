include { DATASETS_SUMMARY }  from '../../../modules/local/datasets_summary'
include { DATASETS_DOWNLOAD } from '../../../modules/local/datasets_download'
include { REFERENCES }        from '../../../modules/local/references'
include { SKANI_SKETCH }      from '../../../modules/local/skanisketch'
include { SKANI_DIST }        from '../../../modules/local/skanidist'
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

            SPESTIMATOR(ch_contigs)

            SPESTIMATOR.out.results
                .map { it -> it [1] }
                .collectFile(
                    storeDir: "${params.outdir}/spestimator/",
                    keepHeader: true,
                    sort: { file -> file.text },
                    name: "spestimator_summary.tsv")
                .set { ch_spestimator_summary }


            ch_versions = ch_versions.mix(SPESTIMATOR.out.versions.first())
            ch_species = ch_species.mix(ch_spestimator_summary )

            SPECIES(ch_species.collect())

            SPECIES.out.species
                .splitText()
                .map{ it -> it.trim()}
                .set{ ch_species_list }

            DATASETS_SUMMARY(ch_species_list.combine(dataset_script))
            ch_versions = ch_versions.mix(DATASETS_SUMMARY.out.versions.first())

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
                    log.info "Successfully identified ${genomeCount} genomes from NCBI for ANI analysis."
                    if (genomeCount == 0) {
                        log.warn "No genomes were identified from NCBI for ANI analysis. This may be due to issues with the NCBI API, or because no reference genomes were identified for the species in the dataset. If you believe there should be reference genomes available, try running this workflow again, or provide references with 'params.reference_genomes'."
                    }
                }

            DATASETS_DOWNLOAD(DATASETS_SUMMARY.out.genomes.mix(SPECIES.out.accessions).collect())

            DATASETS_DOWNLOAD.out.genomes
                .flatten()
                .count()
                .subscribe { genomeCount ->
                    if (genomeCount == 0) {
                        log.warn "No genomes were downloaded from NCBI for ANI analysis. This may be due to issues with the NCBI API, or because no reference genomes were identified for the species in the dataset. If you believe there should be reference genomes available, try running this workflow again, or provide references with 'params.reference_genomes'."
                    } else {
                        log.info "Successfully downloaded ${genomeCount} genomes from NCBI for ANI analysis."
                    }
                }

            ch_reference_genomes = ch_reference_genomes.mix(DATASETS_DOWNLOAD.out.genomes.flatten())
            ch_versions = ch_versions.mix(DATASETS_DOWNLOAD.out.versions)


        } else {
            log.info "Using local files and those provided with this workflow."
            log.info "FYI: To download additional references from NCBI, set 'params.current_datasets'."
            ch_datasets_summary = channel.empty()
        }

        REFERENCES()

        // due to how input files can be added, this removes duplicates
        ch_reference_genomes
            .mix(REFERENCES.out.fastas.flatten())
            .map { file ->
                tuple(file.name, file)
            }
            .groupTuple()
            .map { acc, files ->
                files[0]
            }
            .collect()
            .set { ch_deduplicated_reference_genomes }

        log.info "Running SKANI sketching and distance calculation for average nucleotide identity (ANI) analysis)."
        SKANI_SKETCH(ch_deduplicated_reference_genomes)
        ch_versions = ch_versions.mix(SKANI_SKETCH.out.versions.first())

        SKANI_DIST(ch_contigs, SKANI_SKETCH.out.db)

        SKANI_DIST.out.results
            .map { it -> it [1] }
            .collectFile(
                storeDir: "${params.outdir}/skani/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "skani_summary.tsv")
            .set { ch_skani_summary }

        ch_versions = ch_versions.mix(SKANI_DIST.out.versions.first())

        log.info "Using SKANI results to designate species of contigs"
        SKANI_DIST.out.hits
            .map { meta, contigs, tsv ->
                def lines = tsv.readLines()
                
                if (lines.size() > 1) {
                    def topHitRow = lines[1]
                    def firstCol = topHitRow.split()[0] 
                    def parts = firstCol.split('_')                   
                    return tuple( meta , [parts[0], parts[1]] , contigs )
                } else {
                    // Handle cases where no hits were found
                    return tuple( meta , ["Unknown", "Unknown"] , contigs )
                }
            }
            .view { meta, organism, _ -> 
                "Sample ${meta.id} is predicted to be a ${organism[0]} ${organism[1]}" 
            }
            .set {ch_org_contigs }

        if ( params.msa && ! params.exclude_top_hit ) {
            log.info "Adding top hits from SKANI results to the analysis for multiple sequence alignment (MSA) and phylogenetic analysis. This will add the reference genome with the highest ANI for use in the phylogenetic analysis workflow. If you want to skip this step, set 'params.exclude_top_hit' to true."
            SKANI_DIST.out.top_hit
                .collectFile(name: 'top_hits.txt', newLine: true)
                .splitText()
                .map { it.trim() }
                .filter { it }
                .toList()
                .map { list -> [filterList: list] }
                .set{ ch_skani_top_hits }

            ch_deduplicated_reference_genomes
                .flatten()
                .combine(ch_skani_top_hits)
                .filter { file, meta ->
                    meta.filterList.any { id -> file.name.contains(id) }
                }
                .map { file, meta -> file }
                .view { "Adding SKANI top hit ${it.name} to analysis." }
                .map { it ->
                    def meta = [id:it.baseName]
                    def species = it.name.split("_")[0]
                    def genus = it.name.split("_")[1]
                    tuple( meta, [species, genus], it)
                }
                .set { ch_top_hits }
        } else {
            ch_top_hits = channel.empty()
        }


    emit:
        for_summary      = ch_skani_summary.mix(ch_datasets_summary)
        top_hit          = ch_top_hits
        ch_org_contigs   = ch_org_contigs
        ch_salmonella    = SKANI_DIST.out.salmonella
        ch_ecoli         = SKANI_DIST.out.ecoli
        ch_kleb          = SKANI_DIST.out.kleb
        ch_gas           = SKANI_DIST.out.gas
        ch_strep         = SKANI_DIST.out.strep
        ch_legionella    = SKANI_DIST.out.legionella
        ch_vibrio        = SKANI_DIST.out.vibrio
        ch_acinetobacter = SKANI_DIST.out.acinetobacter
        ch_myco          = SKANI_DIST.out.myco
        ch_gc            = SKANI_DIST.out.gc

        versions    = ch_versions
}

workflow.onComplete {
  log.info "Average nucleotide identity workflow completed at: $workflow.complete"
  log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}