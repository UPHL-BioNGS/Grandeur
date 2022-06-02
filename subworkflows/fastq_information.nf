include { fastqc }                                 from '../modules/fastqc'         addParams(fastq_processes: params.fastq_processes, fastqc_options: params.fastqc_options )
include { mash_sketch; mash_dist }                 from '../modules/mash'           addParams(fastq_processes: params.fastq_processes, mash_reference: params.mash_reference, mash_options: params.mash_options)
include { kraken2_fastq as kraken2 }               from '../modules/kraken2'        addParams(fastq_processes: params.fastq_processes, kraken2_options: params.kraken2_options )
include { lyveset_shuffle; lyveset_cg_pipeline }   from '../modules/lyveset'        addParams(fastq_processes: params.fastq_processes, cg_pipeline_options: params.cg_pipeline_options )
include { serotypefinder_fastq as serotypefinder } from '../modules/serotypefinder' addParams(fastq_processes: params.fastq_processes, serotypefinder_options: params.serotypefinder_options )
include { seqsero2_fastq as seqsero2 }             from '../modules/seqsero2'       addParams(fastq_processes: params.fastq_processes, seqsero2_options: params.seqsero2_options )
include { shigatyper }                             from '../modules/shigatyper'     addParams(fastq_processes: params.fastq_processes, shigatyper_options: params.shigatyper_options )
include { plasmidfinder }                          from '../modules/plasmidfinder'  addParams(fastq_processes: params.fastq_processes, plasmidfinder_options: params.plasmidfinder_options )

workflow fastq_information {
  take:
    reads
    clean_reads
    genome_sizes
    kraken2_db
  main:
    fastqc(reads)
    plasmidfinder(clean_reads)
    mash_sketch(clean_reads)
    mash_dist(mash_sketch.out.files)
    shigatyper(clean_reads.join(mash_dist.out.ecoli_flag, by: 0))
    lyveset_shuffle(clean_reads)

    lyveset_shuffle.out.shuffled
      .join(mash_sketch.out.genome_size, by: 0)
      .join(mash_dist.out.genus, by: 0)
      .join(mash_dist.out.species, by: 0)
      .combine(genome_sizes)
      .set { for_gc }
    lyveset_cg_pipeline(for_gc)

    kraken2(clean_reads.combine(kraken2_db))
    serotypefinder(clean_reads.join(mash_dist.out.ecoli_flag, by: 0))
    seqsero2(clean_reads.join(mash_dist.out.salmonella_flag, by: 0))

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
    shigatyper_predictions  = shigatyper.out.predictions
    shigatyper_cada         = shigatyper.out.cada
    fastqc_1_results        = fastqc.out.fastqc_1_results
    fastqc_2_results        = fastqc.out.fastqc_2_results
    kraken2_top_hit         = kraken2.out.top_hit
    kraken2_top_perc        = kraken2.out.top_perc
    kraken2_top_reads       = kraken2.out.top_reads
    serotypefinder_ogroup   = serotypefinder.out.ogroup
    serotypefinder_hgroup   = serotypefinder.out.hgroup
    seqsero2_profile        = seqsero2.out.profile
    seqsero2_serotype       = seqsero2.out.serotype
    seqsero2_contamination  = seqsero2.out.contamination
    plasmidfinder_hits      = plasmidfinder.out.plasmids
    cg_pipeline_read_length = lyveset_cg_pipeline.out.read_length
    cg_pipeline_quality     = lyveset_cg_pipeline.out.quality
    cg_pipeline_coverage    = lyveset_cg_pipeline.out.coverage
    cg_pipeline_ref_gen_len = lyveset_cg_pipeline.out.ref_genome_length

    // for combined files
    seqsero2_collect        = seqsero2.out.collect

    // for multiqc
    fastqc_multiqc          = fastqc.out.for_multiqc
    kraken2_multiqc         = kraken2.out.for_multiqc
}
