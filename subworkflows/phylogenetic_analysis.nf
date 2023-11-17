include { core_genome_evaluation } from '../modules/grandeur'    addParams(params)
include { heatcluster }            from '../modules/heatcluster' addParams(params)
include { iqtree2 }                from '../modules/iqtree2'     addParams(params)
include { mashtree }               from '../modules/mashtree'    addParams(params)
include { panaroo }                from '../modules/panaroo'     addParams(params)
include { phytreeviz }             from '../modules/phytreeviz'  addParams(params)
include { prokka }                 from '../modules/prokka'      addParams(params)
include { roary }                  from '../modules/roary'       addParams(params)
include { snp_dists }              from '../modules/snp-dists'   addParams(params)

workflow phylogenetic_analysis {
  take:
    evaluat_script
    ch_contigs
    ch_gff
    ch_top_hit

  main:
    for_prokka  = Channel.empty()
    ch_organism = Channel.empty()

    if (params.extras) {
      ch_organism = ch_organism.mix(ch_top_hit.map { it -> if (it) { tuple( it[0] , [ it[1].split("_")[0], it[1].split("_")[1]] )}})

      if ( params.fastani_include ) {
        ch_top_hit
          .map { it -> if (it) { tuple( it[1].split("_", 3)[2], it[2], it[1].split("_")[0, 1]) }}
          .groupTuple(by: 0)
          .map { it -> if (it) { tuple( it[1][0].baseName, it[1][0], it[2][0] ) }}
          .unique()
          .set { ch_representative }

        for_prokka = for_prokka.mix(ch_representative)
      }
    }

    ch_contig_organism = ch_contigs.join( ch_organism, by: 0, remainder: true)

    for_prokka = for_prokka.mix(ch_contig_organism).unique()

    prokka(for_prokka)

    // panaroo(prokka.out.gffs.concat(ch_gff).unique().collect())
    // core_genome_evaluation(panaroo.out.core_gene_alignment.combine(evaluat_script))

    roary(prokka.out.gffs.concat(ch_gff).unique().collect())
    core_genome_evaluation(roary.out.core_gene_alignment.combine(evaluat_script))

    core_genome_evaluation.out.evaluation
      .filter({it[1] as int >= 4})
      .filter({it[2] as int >= params.min_core_genes})
      .map ( it -> it[0] )
      .set{ ch_core_genome }

    iqtree2(ch_core_genome)
    snp_dists(core_genome_evaluation.out.evaluation.map(it->it[0]))
    heatcluster(snp_dists.out.snp_matrix)
    phytreeviz(iqtree2.out.newick)

  emit:
    for_multiqc = prokka.out.for_multiqc.mix(snp_dists.out.snp_matrix).mix(heatcluster.out.for_multiqc).mix(phytreeviz.out.for_multiqc).mix(core_genome_evaluation.out.for_multiqc)
}
