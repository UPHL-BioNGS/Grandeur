//https://github.com/StaPH-B/docker-builds/tree/master/pbptyper/1.0.4
//genotyphi/mykrobe
//shigeifinder

include { mash_dist }                 from '../modules/mash'    addParams(params)
include { mash_sketch }                 from '../modules/mash'  addParams(params)

workflow determine_organism {
  take:
    reads
    clean_reads
    genome_sizes
    fastas

  main:
    mash_sketch(clean_reads)
    mash_dist(mash_sketch.out.files)

    lyveset_shuffle(clean_reads)

    lyveset_shuffle.out.shuffled
      .join(mash_sketch.out.genome_size, by: 0)
      .join(mash_dist.out.genus, by: 0)
      .join(mash_dist.out.species, by: 0)
      .combine(genome_sizes)
      .set { for_gc }
    lyveset_cg_pipeline(for_gc)

    lyveset_cg_pipeline.out.collect
      .collectFile(name: "cg_pipeline_report.txt",
        keepHeader: true,
        sort: true,
        storeDir: "${params.outdir}/cg_pipeline")

  emit:
    // for the summary file
    mash_genome_size        = mash_sketch.out.genome_size
    mash_coverage           = mash_sketch.out.coverage
    mash_species            = mash_dist.out.species
    mash_genus              = mash_dist.out.genus
    mash_full               = mash_dist.out.full
    mash_pvalue             = mash_dist.out.pvalue
    mash_distance           = mash_dist.out.distance
    ecoli_flag              = mash_dist.out.ecoli_flag
    salmonella_flag         = mash_dist.out.salmonella_flag
    klebsiella_flag         = mash_dist.out.klebsiella_flag
    cg_pipeline_read_length = lyveset_cg_pipeline.out.read_length
    cg_pipeline_quality     = lyveset_cg_pipeline.out.quality
    cg_pipeline_coverage    = lyveset_cg_pipeline.out.coverage
    cg_pipeline_ref_gen_len = lyveset_cg_pipeline.out.ref_genome_length

}
