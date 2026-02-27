include { CORE_GENOME_EVALUATION } from '../../../modules/local/core_genome_evaluation'
include { BAKTA }                  from '../../../modules/local/bakta'
include { GOTREE }                 from '../../../modules/local/gotree'
include { HEATCLUSTER }            from '../../../modules/local/heatcluster'
include { IQTREE }                 from '../../../modules/local/iqtree'
include { KSNP4 }                  from '../../../modules/local/ksnp4'
include { MASHTREE }               from '../../../modules/local/mashtree'
include { PANAROO }                from '../../../modules/local/panaroo'
include { PROKKA }                 from '../../../modules/local/prokka'
include { ROARY }                  from '../../../modules/local/roary'
include { SKA2 }                   from '../../../modules/local/ska2'
include { SNPDISTS }               from '../../../modules/local/snpdists'

workflow PHYLOGENETIC_ANALYSIS {
  take:
  evaluat_script
  ch_contigs
  ch_top_hit

  main:
  ch_versions = Channel.empty()
  ch_multiqc  = Channel.empty()
  ch_nwk      = Channel.empty()

  log.info """

Running phylogenetic analysis.

Relevant params and their values:
- 'params.annotator' : ${params.annotator}
    - Specifies what is used to annotate assemblies (impacts what genes are predicted)
    - Options are 'prokka', 'bakta', or 'none'
- 'params.aligner' : ${params.aligner}
    - Specifies what is used to align core genome (the genes that organisms all share)
    - Options are 'panaroo', 'roary', or 'none'
- 'params.min_core_genes' : ${params.min_core_genes}
    - Specifies the minimum number of genes that must be shared.
    - The number of genes that are shared will differ based on the organism and purpose 
      of the analysis
- 'params.min_core_per' : ${params.min_core_per}
    - Specifies the minimum percentage of genes identified that are shared by the other
      input files. If the number of core genes is very low, it may be worth 
      investigating the samples with low shared gene content and considering removing 
      them from the analysis if they are outliers or of poor quality.
    - A general rule of thumb is that at least 70% of the genes should be shared to have 
      a robust core genome when evaluating outbreaks. 
- 'params.exclude_top_hit' : ${params.exclude_top_hit}
    - When 'true', subworkflow does not add the top hits identified in SKANI_DIST
    - When 'false', subworkflow will not add SKANI_DIST results
- 'params.skip_extras' : ${params.skip_extras}
    - When 'true', subworkflow does not run skani and cannot use the skani top hit in 
      analysis

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process           ┃ description                                                        ┃
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ PROKKA            ┃ Predicts genes in assemblies.                                      ┃
┃ BAKTA             ┃ Predicts genes in assemblies.                                      ┃
┃ PANAROO           ┃ Aligns core genomes.                                               ┃
┃ ROARY             ┃ Aligns core genomes.                                               ┃
┃ CORE_GENOME_EV... ┃ Custom process to identify number of core genes identified and the ┃
┃                   ┃ percentage of how many identified genes are in this core for each  ┃
┃                   ┃ input.                                                             ┃
┃ KSNP4             ┃ Uses a k-mer based approach to identify core SNPs.                 ┃
┃ MASHTREE          ┃ Uses a k-mer based approach to calculate distances between genomes ┃
┃                   ┃ and construct a tree.                                              ┃
┃ SKA2              ┃ Uses a k-mer based approach to align core genes.                   ┃
┃ IQTREE            ┃ Uses a maximum likelihood approach to construct a tree from a core ┃
┃                   ┃ outputs from SKA2 and core genome alignment.                       ┃
┃ GOTREE            ┃ Creating visualizations of the trees with GOTREE.                  ┃
┃ SNPDISTS          ┃ Calculating SNP distance matrix with SNPDISTS from gene alignments.┃
┃ HEATCLUSTER       ┃ Creates a clustered heatmap from SNP distance matrices.            ┃
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

  // adding in organism and top ani hit
  if ( ! params.skip_extras ) {
    ch_organism = ch_top_hit.map { it -> if (it) { tuple( it[0] , [ it[1].split("_")[0], it[1].split("_")[1]] )}}

    if ( ! params.exclude_top_hit ) {
      ch_top_hit
        .map { it -> if (it) { tuple( it[1].split("_", 3)[2], it[2], it[1].split("_")[0, 1]) }}
        .groupTuple(by: 0)
        .map { it -> 
          if (it) {
            def meta = [id:it[1][0].baseName] 
            tuple( meta, it[1][0], it[2][0] ) }}
        .unique()
        .set { ch_representative }

      ch_preannotation = ch_contigs.join( ch_organism, by: 0, remainder: true).mix(ch_representative).filter{it}
    } else {
      ch_preannotation = ch_contigs.join( ch_organism, by: 0, remainder: true)
    }
  } else {
    // skipping ani and top hit
    ch_preannotation  = ch_contigs.map{ it -> tuple(it[0], it[1], null)}
  }
  
  if (params.annotator == 'prokka' ) {
    PROKKA(ch_preannotation.unique())
    
    ch_versions = ch_versions.mix(PROKKA.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(PROKKA.out.for_multiqc)
    ch_gff      = PROKKA.out.gff
  } else if (params.annotator == 'bakta') {
    BAKTA(ch_preannotation.unique())
    
    ch_versions = ch_versions.mix(BAKTA.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(BAKTA.out.for_multiqc)
    ch_gff      = BAKTA.out.gff

  } else {
    ch_gff = Channel.empty()

  }

  if (params.aligner == 'panaroo') {
    PANAROO(ch_gff.unique().collect())

    ch_core     = PANAROO.out.core_gene_alignment
    ch_versions = ch_versions.mix(PANAROO.out.versions)

  } else if (params.aligner == 'roary') {
    ROARY(ch_gff.unique().collect())

    ch_core     = ROARY.out.core_gene_alignment
    ch_versions = ch_versions.mix(ROARY.out.versions)
  } else {
    ch_core     = Channel.empty()
  }

  CORE_GENOME_EVALUATION(ch_core.combine(evaluat_script))

  CORE_GENOME_EVALUATION.out.evaluation
    .splitText()
    .first()
    .map { it.trim().split(',') }
    .view { num_samples, num_core_genes, core_genome_per ->
        "Core Genome Evaluation Complete: Found ${num_core_genes} core genes (Core Percentage: ${core_genome_per}%)"
    }
    .tap { ch_stats ->
        ch_stats.subscribe { num_samples, num_core_genes, core_genome_per ->
            if (params.min_core_genes && (num_core_genes as int) < params.min_core_genes) {
                log.warn "WARNING: Core genes (${num_core_genes}) is below the minimum threshold of ${params.min_core_genes}!"
            }
            if (params.min_core_per && (core_genome_per as float) < params.min_core_per) {
                log.warn "WARNING: Core percentage (${core_genome_per}%) is below the minimum threshold of ${params.min_core_per}%!"
            }
        }
    }
    .combine(ch_core)
    // filter out if there are too few core genes or a very low core genome percentage.
    .filter { num_samples, num_core_genes, core_genome_per, core_files ->
        def pass_genes = params.min_core_genes ? (num_core_genes as int >= params.min_core_genes) : true
        def pass_per   = params.min_core_per   ? (core_genome_per as float >= params.min_core_per) : true
        return pass_genes && pass_per
    }
    .map{ it -> it[-2]}
    .set { ch_core_genome }

  ch_multiqc = ch_multiqc.mix(CORE_GENOME_EVALUATION.out.for_multiqc)

  KSNP4(ch_contigs.combine(ch_top_hit))
  ch_nwk = ch_nwk.mix(KSNP4.out.newick)
  ch_versions = ch_versions.mix(KSNP4.out.versions.first())

  MASHTREE(ch_preannotation.map{it -> if (it) { tuple( it[1]) }}.collect())
  ch_nwk = ch_nwk.mix(MASHTREE.out.newick)
  ch_versions = ch_versions.mix(MASHTREE.out.versions)

  SKA2(ch_contigs.combine(ch_top_hit))
  ch_versions = ch_versions.mix(SKA2.out.versions.first())
    
  IQTREE(ch_core_genome.mix(SKA2.out.aln))
  ch_nwk = ch_nwk.mix(IQTREE.out.newick)
  ch_versions = ch_versions.mix(IQTREE.out.versions.first())

  GOTREE(ch_nwk)
  ch_versions = ch_versions.mix(GOTREE.out.versions.first())
  ch_multiqc  = ch_multiqc.mix(GOTREE.out.for_multiqc)

  SNPDISTS(ch_core_genome.mix(ska2.out.aln))
  ch_versions = ch_versions.mix(SNPDISTS.out.versions)
  ch_multiqc  = ch_multiqc.mix(SNPDISTS.out.snp_matrix)

  HEATCLUSTER(SNPDISTS.out.snp_matrix)
  ch_versions = ch_versions.mix(HEATCLUSTER.out.versions)
  ch_multiqc  = ch_multiqc.mix(HEATCLUSTER.out.for_multiqc)

  emit:
  for_multiqc = ch_multiqc
  versions    = ch_versions
}

if ( params.msa ) {
    workflow.onComplete {
        log.info """------------------------------------------------------

PHYLOGENETIC ANALYSIS subworkflow completed at: $workflow.complete

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Subworkflow Output Files                              ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│   'params.outdir'                                     │"""
        if ( params.aligner ) {
    log.info """│    ├── spestimator                                    │
│    │   └── spestimator_summary.tsv                    │
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
