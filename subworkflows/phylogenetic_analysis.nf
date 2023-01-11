include { iqtree2 }        from '../modules/iqtree2'   addParams(params)
include { prokka }         from '../modules/prokka'    addParams(params)
include { representative } from '../modules/grandeur'  addParams(params)
include { roary }          from '../modules/roary'     addParams(params)
include { snp_dists }      from '../modules/snp-dists' addParams(params)

workflow phylogenetic_analysis {
  take:
    ch_contigs
    ch_gff
    ch_organism
    ch_genomes

  main:
    representative(ch_organism.map{ it -> it[1][2]}.unique().combine(ch_genomes))

    ch_representative = representative.out.representative.map{ it -> tuple(it[3], it[0] , [it[1], it[2], it[3]] )}

    ch_contigs
      .join(ch_organism, by: 0, remainder: true)
      .mix(ch_representative)
      .set { for_prokka }

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
