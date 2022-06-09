include { roary }     from '../modules/roary'     addParams(phylogenetic_processes: params.phylogenetic_processes, roary_options: params.roary_options )
include { iqtree2 }   from '../modules/iqtree2'   addParams(phylogenetic_processes: params.phylogenetic_processes, iqtree2_options: params.iqtree2_options )
include { snp_dists } from '../modules/snp-dists' addParams(phylogenetic_processes: params.phylogenetic_processes, snp_dists_options: params.snp_dists_options )
include { prokka }    from '../modules/prokka'    addParams(phylogenetic_processes: params.phylogenetic_processes, prokka_options: params.prokka_options, outgroup: params.outgroup )

workflow phylogenetic_analysis {
  take:
    contigs
    mash_species
    mash_genus
    gffs
  main:
    contigs
      .join(mash_genus, by: 0, remainder: true)
      .join(mash_species, by: 0, remainder: true)
      .set { for_prokka }
    prokka(for_prokka)
    roary(prokka.out.gffs.concat(gffs).collect())
    iqtree2(roary.out.core_gene_alignment)
    snp_dists(roary.out.core_gene_alignment)
  emit:
    prokka_multiqc = prokka.out.for_multiqc
}
