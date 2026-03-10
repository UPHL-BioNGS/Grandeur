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
        log.info """

Running average nucleotide identity (ANI) analysis).

More information can be found at https://github.com/UPHL-BioNGS/Grandeur/wiki/average_nucleotide_identity.

Relevant params and their values:
- 'params.current_datasets' : ${params.current_datasets}
    - When 'true', subworkflow will download additional references from NCBI
    - When 'false', subworkflow will use local references and SPESTIMATOR, SPECIES,
      DATASETS_SUMMARY, and DATASETS_DOWNLOAD will be skipped.
    - Downloading reference genomes uses a third-party API that is not controlled by 
      the Grandeur developers, requires the workflow to have internet access, and may be 
      slow or have issues.

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process           ┃ description                                                        ┃
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ SPESTIMATOR       ┃ Uses 16S to identify range of species to download.                 ┃
┃ SPECIES           ┃ Reference genomes are identified using results from SPESTIMATOR    ┃
┃                   ┃ and MASH, as well as SYLPH and KRAKEN2 (if their respective        ┃
┃                   ┃ databases are provided).                                           ┃
┃ DATASETS_SUMMARY  ┃ Looks up refence accession for each identified species.            ┃
┃ DATASETS_DOWNLOAD ┃ Downloads reference genomes.                                       ┃
┃ REFERENCES        ┃ Loads stored reference genomes.                                    ┃
┃ SKANI_SKETCH      ┃ Creates a sketch of all reference genomes.                         ┃
┃ SKANI_DIST        ┃ Estimates distance of input files to references. Core process of   ┃
┃                   ┃ organism estimation.                                               ┃ 
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

        ch_versions = channel.empty()
        ch_summary = channel.empty()

        if ( params.current_datasets ) {

            SPESTIMATOR(ch_contigs)

            SPESTIMATOR.out.results
                .map { it -> it [1] }
                .collectFile(
                    storeDir: "${params.outdir}/spestimator/",
                    keepHeader: true,
                    sort: { file -> file.text },
                    name: "spestimator_summary.csv")
                .set { ch_spestimator_summary }


            ch_versions = ch_versions.mix(SPESTIMATOR.out.versions.first())
            ch_summary  = ch_summary.mix(ch_spestimator_summary)
            ch_species = ch_species.mix(ch_spestimator_summary )

            // could be a channel, but some mash results are very long and may overload headnodes
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

            ch_summary = ch_summary.mix(ch_datasets_summary)

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

        SKANI_SKETCH(ch_deduplicated_reference_genomes)
        ch_versions = ch_versions.mix(SKANI_SKETCH.out.versions)

        SKANI_DIST(ch_contigs, SKANI_SKETCH.out.db)

        SKANI_DIST.out.skani
            .collectFile(
                storeDir: "${params.outdir}/skani/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "skani_summary.tsv")
            .set { ch_skani_summary }

        ch_summary  = ch_summary.mix(ch_skani_summary)
        ch_versions = ch_versions.mix(SKANI_DIST.out.versions.first())

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
            log.info "Adding top hits from SKANI results to the analysis."
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
        for_summary      = ch_summary
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

        versions         = ch_versions
}

if ( ! params.skip_extras ) {
    workflow.onComplete {
        log.info """------------------------------------------------------

AVERAGE NUCLEOTIDE IDENTITY subworkflow completed at: $workflow.complete

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Subworkflow Output Files                              ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│   'params.outdir'                                     │"""
        if ( params.current_datasets ) {
    log.info """│    ├── spestimator                                    │
│    │   └── spestimator_summary.csv                    │
│    ├── datasets                                       │
│    │   └── datasets_summary.csv                       │"""
        }

        log.info """│    └── skani                                          │
│        └── skani_summary.tsv                          │
└───────────────────────────────────────────────────────┘

------------------------------------------------------
"""
    }
}
