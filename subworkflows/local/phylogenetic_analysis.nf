include { CORE_GENOME_EVALUATION } from '../../modules/local/local'
include { BAKTA }                  from '../../modules/local/bakta'
include { HEATCLUSTER }            from '../../modules/local/heatcluster'
include { IQTREE2 }                from '../../modules/local/iqtree2'
include { MASHTREE }               from '../../modules/local/mashtree'
include { PANAROO }                from '../../modules/local/panaroo'
include { PHYTREEVIZ }             from '../../modules/local/phytreeviz'
include { PROKKA }                 from '../../modules/local/prokka'
include { ROARY }                  from '../../modules/local/roary'
include { SNPDISTS }               from '../../modules/local/snp-dists'

workflow PHYLOGENETIC_ANALYSIS {
  take:
  evaluat_script
  ch_contigs
  ch_top_hit

  main:
  ch_versions = Channel.empty()
  ch_multiqc  = Channel.empty()

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

      ch_preannotation = ch_contigs.join( ch_organism, by: 0, remainder: true).mix(ch_representative)
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
    .map{ it -> it.trim()}
    .map { it ->
      def (num_samples, num_core_genes, core_genome_per) = it.split(',')
      return [num_samples, num_core_genes, core_genome_per]
    }
    .combine(ch_core)
    .set { ch_core_genome }

  if (params.min_core_genes) {
    ch_core_genome = ch_core_genome.filter{it[1] as int >= params.min_core_genes}
  }

  if (params.min_core_per) {
    ch_core_genome = ch_core_genome.filter{it[2] as float >= params.min_core_per}
  }

  ch_core_genome = ch_core_genome.map{ it -> it[-2]}

  ch_multiqc = ch_multiqc.mix(CORE_GENOME_EVALUATION.out.for_multiqc)

  // TODO : if channel doesn't go to to iqtree2, then send to mashtree

  // phylogenetic trees
  MASHTREE(ch_preannotation.map{it -> if (it) { tuple( it[1]) }}.collect())
  ch_versions = ch_versions.mix(MASHTREE.out.versions)
    
  IQTREE2(ch_core_genome)
  ch_versions = ch_versions.mix(IQTREE2.out.versions)

  PHYTREEVIZ(IQTREE2.out.newick.mix(MASHTREE.out.newick))
  ch_versions = ch_versions.mix(PHYTREEVIZ.out.versions.first())
  ch_multiqc  = ch_multiqc.mix(PHYTREEVIZ.out.for_multiqc)

  // SNP matrix
  SNPDISTS(ch_core_genome)
  ch_versions = ch_versions.mix(SNPDISTS.out.versions)
  ch_multiqc  = ch_multiqc.mix(SNPDISTS.out.snp_matrix)

  HEATCLUSTER(SNPDISTS.out.snp_matrix)
  ch_versions = ch_versions.mix(HEATCLUSTER.out.versions)
  ch_multiqc  = ch_multiqc.mix(HEATCLUSTER.out.for_multiqc)

  emit:
  for_multiqc = ch_multiqc
  versions    = ch_versions
}
