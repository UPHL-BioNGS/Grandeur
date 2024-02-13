include { core_genome_evaluation } from '../modules/local/local'       addParams(params)
include { heatcluster }            from '../modules/local/heatcluster' addParams(params)
include { iqtree2 }                from '../modules/local/iqtree2'     addParams(params)
include { mashtree }               from '../modules/local/mashtree'    addParams(params)
include { panaroo }                from '../modules/local/panaroo'     addParams(params)
include { phytreeviz }             from '../modules/local/phytreeviz'  addParams(params)
include { prokka }                 from '../modules/local/prokka'      addParams(params)
include { roary }                  from '../modules/local/roary'       addParams(params)
include { snp_dists }              from '../modules/local/snp-dists'   addParams(params)

workflow phylogenetic_analysis {
  take:
    evaluat_script
    ch_contigs
    ch_top_hit

  main:

    // adding in organism and top ani hit
    if ( ! params.skip_extras ) {
      ch_organism = ch_top_hit.map { it -> if (it) { tuple( it[0] , [ it[1].split("_")[0], it[1].split("_")[1]] )}}

      if ( ! params.exclude_top_hit ) {
        ch_top_hit
          .map { it -> if (it) { tuple( it[1].split("_", 3)[2], it[2], it[1].split("_")[0, 1]) }}
          .groupTuple(by: 0)
          .map { it -> 
            if (it) {
              meta = [id:it[1][0].baseName] 
              tuple( meta, it[1][0], it[2][0] ) }}
          .unique()
          .set { ch_representative }

        for_prokka = ch_contigs.join( ch_organism, by: 0, remainder: true).mix(ch_representative)
      } else {
        for_prokka = ch_contigs.join( ch_organism, by: 0, remainder: true)
      }
    } else {
      // skipping ani and top hit
      ch_organism = Channel.empty() 
      for_prokka  = ch_contigs.map{ it -> tuple(it[0], it[1], null)}
    }

    prokka(for_prokka.unique())

    panaroo(prokka.out.gff.unique().collect())
    core_genome_evaluation(panaroo.out.core_gene_alignment.combine(evaluat_script))

    core_genome_evaluation.out.evaluation
      .filter({it[1] as int >= 4})
      .filter({it[2] as int >= params.min_core_genes})
      .map ( it -> it[0] )
      .set{ ch_core_genome }

    // TODO : if channel doesn't go to to iqtree2, then send to mashtree

    // phylogenetic trees
    mashtree(for_prokka.map{it -> tuple( it[1]) }.collect())
    iqtree2(ch_core_genome)
    phytreeviz(iqtree2.out.newick.mix(mashtree.out.newick))

    // SNP matrix
    snp_dists(core_genome_evaluation.out.evaluation.map(it->it[0]))
    heatcluster(snp_dists.out.snp_matrix)

  emit:
    for_multiqc = prokka.out.for_multiqc.mix(snp_dists.out.snp_matrix).mix(heatcluster.out.for_multiqc).mix(phytreeviz.out.for_multiqc).mix(core_genome_evaluation.out.for_multiqc)
    versions    = prokka.out.versions.first().mix(panaroo.out.versions).mix(mashtree.out.versions).mix(iqtree2.out.versions).mix(phytreeviz.out.versions.first()).mix(snp_dists.out.versions).mix(heatcluster.out.versions)
}
