include { amrfinderplus }                           from '../modules/amrfinderplus'  addParams(contig_processes: params.contig_processes, amrfinderplus_options: params.amrfinderplus_options )
include { fastani }                                 from '../modules/fastani'        addParams(contig_processes: params.contig_processes, fastani_options: params.fastani_options )
include { kleborate }                               from '../modules/kleborate'      addParams(contig_processes: params.contig_processes, kleborate_options: params.kleborate_options )
include { kraken2_fasta as kraken2 }                from '../modules/kraken2'        addParams(contig_processes: params.contig_processes, kraken2_options: params.kraken2_options )
include { mlst }                                    from '../modules/mlst'           addParams(contig_processes: params.contig_processes, mlst_options: params.mlst_options )
include { quast }                                   from '../modules/quast'          addParams(contig_processes: params.contig_processes, quast_options: params.quast_options )
include { seqsero2_fasta as seqsero2 }              from '../modules/seqsero2'       addParams(contig_processes: params.contig_processes, seqsero2_options: params.seqsero2_options )
include { serotypefinder_fasta as serotypefinder }  from '../modules/serotypefinder' addParams(contig_processes: params.contig_processes, serotypefinder_options: params.serotypefinder_options )
include { plasmidfinder }                           from '../modules/plasmidfinder'  addParams(contig_processes: params.contig_processes, plasmidfinder_options: params.plasmidfinder_options )

workflow contig_information {
  take:
    contigs
    mash_species
    mash_genus
    salmonella_flag
    ecoli_flag
    klebsiella_flag
    fastani_genomes
    kraken2_db
  main:
    mlst(contigs)
    quast(contigs)
    plasmidfinder(contigs)
    fastani(contigs.combine(fastani_genomes))
    kleborate(contigs.join(klebsiella_flag, by:0))
    seqsero2(contigs.join(salmonella_flag,  by:0))
    serotypefinder(contigs.join(ecoli_flag, by:0))
    kraken2(contigs.combine(kraken2_db))

    contigs
      .join(mash_genus,                     by: 0)
      .join(mash_species,                   by: 0)
      .set {for_amrfinder}
    amrfinderplus(for_amrfinder)

    mlst.out.collect
      .collectFile(name: "mlst_result.tsv",
        keepHeader: false,
        sort: true,
        storeDir: "${params.outdir}/mlst")

    kleborate.out.collect
      .collectFile(name: "kleborate_results.txt",
        keepHeader: true,
        sort: true,
        storeDir: "${params.outdir}/kleborate")

    amrfinderplus.out.collect
      .collectFile(name: "ncbi-AMRFinderplus.txt",
        keepHeader: true,
        sort: true,
        storeDir: "${params.outdir}/ncbi-AMRFinderplus")

    fastani.out.collect
      .collectFile(name: "fastani.out.txt",
        sort: true,
        storeDir: "${params.outdir}/fastani")

    quast.out.collect
      .collectFile(name: "report.tsv",
        keepHeader: true,
        sort: true,
        storeDir: "${params.outdir}/quast")

  emit:
    // for summary
    amrfinder_amr_genes     = amrfinderplus.out.amr_genes
    amrfinder_vir_genes     = amrfinderplus.out.vir_genes
    fastani_ref             = fastani.out.ref
    fastani_ani_score       = fastani.out.ani
    fastani_fragment        = fastani.out.fragment
    fastani_total           = fastani.out.total
    kleborate_score         = kleborate.out.score
    kleborate_mlst          = kleborate.out.mlst
    mlst_sttype             = mlst.out.mlst
    kraken2_top_hit         = kraken2.out.top_hit
    kraken2_top_perc        = kraken2.out.top_perc
    kraken2_top_reads       = kraken2.out.top_reads
    serotypefinder_ogroup   = serotypefinder.out.ogroup
    serotypefinder_hgroup   = serotypefinder.out.hgroup
    seqsero2_profile        = seqsero2.out.profile
    seqsero2_serotype       = seqsero2.out.serotype
    seqsero2_contamination  = seqsero2.out.contamination
    quast_gc                = quast.out.gc
    quast_contigs           = quast.out.contigs
    quast_nfifty            = quast.out.nfifty
    quast_length            = quast.out.length
    plasmidfinder_hits      = plasmidfinder.out.plasmids

    // for combined files
    seqsero2_collect        = seqsero2.out.collect

    // for multiqc
    kraken2_multiqc         = kraken2.out.for_multiqc
    quast_multiqc           = quast.out.for_multiqc
}
