include { iqtree2 }        from '../modules/iqtree2'   addParams(params)
include { prokka }         from '../modules/prokka'    addParams(params)
include { representative } from '../modules/grandeur'  addParams(params)
include { roary }          from '../modules/roary'     addParams(params)
include { snp_dists }      from '../modules/snp-dists' addParams(params)

workflow phylogenetic_analysis {
  take:
    ch_contigs
    ch_gff
    ch_top_hit

  main:
    for_prokka = Channel.empty()
    if (params.extras) {
      ch_top_hit
        .map { it -> tuple( it[0] , [ it[1].split("_")[0], it[1].split("_")[1]] )}
        .set { ch_organism }

      if ( params.fastani_include ) {
        ch_top_hit
          .map { it -> tuple( it[1].split("_", 3)[2].replaceAll(~/.fna/,""), it[2], it[1].split("_")[0, 1]) }
          .groupTuple(by: 0)
          .map { it -> tuple( it[0], it[1][0], it[2][0] ) }
          .unique()
          .set { ch_representative }

        for_prokka = for_prokka.mix(ch_representative)
      }
    } else {
      ch_organism = Channel.empty()
    }

    ch_contig_organism = ch_contigs.join( ch_organism, by: 0, remainder: true)

    for_prokka = for_prokka.mix(ch_contig_organism)

    prokka( for_prokka )
    roary(prokka.out.gffs.concat(ch_gff).collect())

    roary.out.core_gene_alignment
      .filter ({ it[1] as int >= 4 })
      .filter ({ it[2] as int >= params.roary_min_genes })
      .map ( it -> it[0] )
      .set{ ch_core_genome }

    iqtree2(ch_core_genome)
    snp_dists(roary.out.core_gene_alignment.map( it -> it[0] ))

  emit:
    for_multiqc = prokka.out.for_multiqc
}
