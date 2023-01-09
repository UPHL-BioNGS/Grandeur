include { roary }     from '../modules/roary'     addParams(params)
include { iqtree2 }   from '../modules/iqtree2'   addParams(params)
include { snp_dists } from '../modules/snp-dists' addParams(params)
include { prokka }    from '../modules/prokka'    addParams(params)

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
