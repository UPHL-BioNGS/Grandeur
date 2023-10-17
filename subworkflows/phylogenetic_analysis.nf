include { iqtree2 }            from '../modules/iqtree2'   addParams(params)
include { prokka }             from '../modules/prokka'    addParams(params)
include { roary }              from '../modules/roary'     addParams(params)
include { snp_dists }          from '../modules/snp-dists' addParams(params)
include { snp_matrix_heatmap } from '../modules/grandeur'  addParams(params)

workflow phylogenetic_analysis {
  take:
    snpmtrx_script
    ch_contigs
    ch_gff
    ch_top_hit

  main:
    for_prokka  = Channel.empty()
    ch_organism = Channel.empty()

    if (params.extras) {
      ch_organism = ch_organism.mix(ch_top_hit.map { it -> tuple( it[0] , [ it[1].split("_")[0], it[1].split("_")[1]] )})

      if ( params.fastani_include ) {
        ch_top_hit
          .map { it -> if (it) { tuple( it[1].split("_", 3)[2], it[2], it[1].split("_")[0, 1]) }}
          .groupTuple(by: 0)
          .map { it -> tuple( it[1][0].baseName, it[1][0], it[2][0] ) }
          .unique()
          .set { ch_representative }

        for_prokka = for_prokka.mix(ch_representative)
      }
    }

    ch_contig_organism = ch_contigs.join( ch_organism, by: 0, remainder: true)

    for_prokka = for_prokka.mix(ch_contig_organism).unique()

    prokka(for_prokka)

    roary(prokka.out.gffs.concat(ch_gff).unique().collect())

    roary.out.core_gene_alignment
      .filter ({ it[1] as int >= 4 })
      .filter ({ it[2] as int >= params.roary_min_genes })
      .map ( it -> it[0] )
      .set{ ch_core_genome }

    iqtree2(ch_core_genome)
    snp_dists(roary.out.core_gene_alignment)
    snp_matrix_heatmap(snp_dists.out.snp_matrix.combine(snpmtrx_script))

  emit:
    for_multiqc = prokka.out.for_multiqc.mix(snp_dists.out.for_multiqc).mix(snp_matrix_heatmap.out.for_multiqc)
}
