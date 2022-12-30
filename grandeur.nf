#!/usr/bin/env nextflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Welcome to this workflow! Issues and contributions are gladly accepted at https://github.com/UPHL-BioNGS/Grandeur .

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

println("Currently using the Grandeur workflow for use with microbial sequencing. The view is great from 8299 feet (2530 meters) above sea level.\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: ${workflow.manifest.version}")
println("")

nextflow.enable.dsl               = 2

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Defining params

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

params.outdir                     = 'grandeur'
params.maxcpus                    = 12
params.medcpus                    = 4
params.minimum_reads              = 10000

// getting some tests going
params.test                       = false

// connecting to phoenix
params.phoenix_wf                 = false

// connecting to Donut Falls
params.donut_falls_wf             = false

// input files
params.reads                      = workflow.launchDir + '/reads'
params.fastas                     = workflow.launchDir + '/fastas'
params.gff                        = workflow.launchDir + '/gff'

// external files
params.kraken2_db                 = ''
params.blast_db                   = ''
params.mash_ref                   = '/db/RefSeqSketchesDefaults.msh'
params.plasmidfinder_ref          = ''
params.fastani_ref                = workflow.projectDir + "/configs/fastani_ref.tar.gz"
params.genome_sizes               = workflow.projectDir + "/configs/genome_sizes.json"

// tool-specific command line options
params.amrfinderplus_options      = ''
params.bbduk_options              = 'k=31 hdist=1'
params.bbmap_options              = ''
params.blast_db_type              = 'nt'
params.blastn_options             = '-max_target_seqs 10 -max_hsps 1 -evalue 1e-25'
params.blobtools_create_options   = ''
params.blobtools_view_options     = ''
params.blobtools_plot_options     = '--format png -r species'
params.blobtools_bbmap_options    = ''
params.cg_pipeline_options        = '--qual_offset 33 --minLength 1'
params.current_datasets           = true
params.datasets_max_genomes       = 5
params.fastani_options            = '--matrix'
params.fastp_options              = '--detect_adapter_for_pe'
params.fastqc_options             = ''
params.iqtree2_options            = '-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000'
params.iqtree2_outgroup           = ''
params.kleborate_options          = '-all'
params.kraken2_options            = ''
params.mash_sketch_options        = '-m 2'
params.mash_dist_options          = '-v 0 -d 0.5'
params.mlst_options               = ''
params.multiqc_options            = ''
params.plasmidfinder_options      = ''
params.prokka_options             = '--mincontiglen 500 --compliant --locustag locus_tag --centre STAPHB'
params.quast_options              = ''
params.roary_options              = ''
params.seqsero2_options           = '-m a -b mem'
params.serotypefinder_options     = ''
params.snp_dists_options          = '-c'
params.spades_options             = '--isolate'

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Sharing params with subworkflows

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

include { de_novo_alignment }           from './subworkflows/de_novo_alignment.nf'            addParams(params)
include { blobtools }                   from './subworkflows/blobtools.nf'                    addParams(params)
include { kraken2 }                     from './subworkflows/kraken2.nf'                      addParams(params)
include { average_nucleotide_identity } from './subworkflows/average_nucleotide_identity.nf'  addParams(params)

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for input files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Getting the fastq files
Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                  "${params.reads}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .view { "Paired-end fastq files found : ${it[0]}" }
  .unique()
  .set { ch_reads }

// Getting contig or fasta files
Channel
  .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
  .map { file -> tuple(file.baseName, file) }
  .unique()
  .set { ch_fastas }

// Getting fasta files that have been annotated with prokka
Channel.fromPath("${params.gff}/*.gff", type: 'file')
  .view { "gff file : $it" }
  .unique()
  .set { ch_gffs }

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for database files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


// Getting the file with genome sizes of common organisms for cg-pipeline. The End User can use their own file and set with a param
ch_genome_sizes = Channel.fromPath(params.genome_sizes, type:'file')

// Getting the reference genomes for fastANI
ch_fastani_genomes = Channel.fromPath("${params.fastani_ref}", type:'file')

// Getting the database for blobtools
ch_blast_db = params.blast_db
              ? Channel
                  .fromPath(params.blast_db, type:'dir')
                  .ifEmpty{
                    println("No blast database was found at ${params.blast_db}")
                    println("Set 'params.blast_db' to directory with blast database")
                    exit 1
                  }
                  .view { "Local Blast Database for Blobtools : $it" }
              : Channel.empty()

// Getting the kraken2 database
ch_kraken2_db = params.kraken2_db
              ? Channel
                  .fromPath(params.kraken2_db, type:'dir')
                  .ifEmpty{
                    println("No kraken2 database was found at ${params.kraken2_db}")
                    println("Set 'params.kraken2_db' to directory with kraken2 database")
                    exit 1
                  }
                  .view { "Local kraken2 database : $it" }
              : Channel.empty()

println("The files and directory for results is " + params.outdir)
println("The maximum number of CPUS for any one process is ${params.maxcpus}")

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Workflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


workflow {

  ch_for_multiqc = Channel.empty()

  if ( params.test )            { get_test_files() }
  if ( params.phoenix_wf )      { phoenix(input) }
  if ( params.donut_falls_wf )  { donut_falls(input) }

  // either phoenix or de_novo_alignmnet is required
  de_novo_alignment(ch_reads)

  CONTIGS = de_novo_alignment.out.contigs
  READS   = de_novo_alignment.out.clean_reads

  // optional subworkflow blobtools (useful for interspecies contamination)
  if (params.blast_db) {
    blobtools(READS, CONTIGS, ch_blast_db)
    blobtools_results = blobtools.out.for_summary
  } else {
    blobtools_results = Channel.empty()
  }

  // optional subworkflow kraken2 (useful for interspecies contamination)
  if (params.kraken2_db) {
    kraken2(READS, CONTIGS, ch_kraken2_db)
    ch_for_multiqc    = ch_for_multiqc.mix(kraken2.out.for_multiqc)
    kraken2_results   = kraken2.out.for_summary
  } else {
    kraken2_results   = Channel.empty()
  }

  // determining organisms in sample
  average_nucleotide_identity(
    kraken2_results
    .concat(blobtools_results)
    .collect(),
    CONTIGS,
    ch_fastani_genomes
  ) 

  // optional subworkflow for comparing shared genes




  // creating txt summary file and multiqc report

  // contigs             = de_novo_alignment.out.contigs.mix(fastas)
  // mash(fastas)
  // mash_species        = fastq_information.out.mash_species.mix(mash.out.species)
  // mash_genus          = fastq_information.out.mash_genus.mix(mash.out.genus)
  // salmonella_flag     = fastq_information.out.salmonella_flag.mix(mash.out.salmonella_flag)
  // ecoli_flag          = fastq_information.out.ecoli_flag.mix(mash.out.ecoli_flag)
  // klebsiella_flag     = fastq_information.out.klebsiella_flag.mix(mash.out.klebsiella_flag)

  // contig_information(contigs, mash_species, mash_genus, salmonella_flag, ecoli_flag, klebsiella_flag, fastani_genomes, local_kraken2)

  // phylogenetic_analysis(contigs, mash_species, mash_genus, gffs )

  // multiqc(de_novo_alignment.out.fastp_multiqc.collect().ifEmpty([]),
  //         de_novo_alignment.out.bbduk_multiqc.collect().ifEmpty([]),
  //         contig_information.out.kraken2_multiqc.collect().ifEmpty([]),
  //         contig_information.out.quast_multiqc.collect().ifEmpty([]),
  //         fastq_information.out.fastqc_multiqc.collect().ifEmpty([]),
  //         fastq_information.out.kraken2_multiqc.collect().ifEmpty([]),
  //         phylogenetic_analysis.out.prokka_multiqc.collect().ifEmpty([]))

  // reads
  //   .mix(fastas)
  //   // de_novo_alignment
  //   .join(de_novo_alignment.out.phix_reads                                                                , remainder: true, by: 0)
  //   .join(de_novo_alignment.out.fastp_results                                                             , remainder: true, by: 0)

  //   // fastq_information
  //   .join(fastq_information.out.fastqc_1_results                                                          , remainder: true, by: 0)
  //   .join(fastq_information.out.fastqc_2_results                                                          , remainder: true, by: 0)
  //   .join(fastq_information.out.cg_pipeline_read_length                                                   , remainder: true, by: 0)
  //   .join(fastq_information.out.cg_pipeline_quality                                                       , remainder: true, by: 0)
  //   .join(fastq_information.out.cg_pipeline_coverage                                                      , remainder: true, by: 0)
  //   .join(fastq_information.out.cg_pipeline_ref_gen_len                                                   , remainder: true, by: 0)
  //   .join(fastq_information.out.shigatyper_predictions                                                    , remainder: true, by: 0)
  //   .join(fastq_information.out.shigatyper_cada                                                           , remainder: true, by: 0)
  //   .join(fastq_information.out.kraken2_top_hit                                                           , remainder: true, by: 0)
  //   .join(fastq_information.out.kraken2_top_perc                                                          , remainder: true, by: 0)
  //   .join(fastq_information.out.kraken2_top_reads                                                         , remainder: true, by: 0)

  //   // mash
  //   .join(fastq_information.out.mash_genome_size                                                          , remainder: true, by: 0)
  //   .join(fastq_information.out.mash_coverage                                                             , remainder: true, by: 0)
  //   .join(mash_genus                                                                                      , remainder: true, by: 0)
  //   .join(mash_species                                                                                    , remainder: true, by: 0)
  //   .join(fastq_information.out.mash_full.mix(mash.out.full)                                              , remainder: true, by: 0)
  //   .join(fastq_information.out.mash_pvalue.mix(mash.out.pvalue)                                          , remainder: true, by: 0)
  //   .join(fastq_information.out.mash_distance.mix(mash.out.distance)                                      , remainder: true, by: 0)

  //   // contig_information
  //   .join(contig_information.out.seqsero2_profile                                                         , remainder: true, by: 0)
  //   .join(contig_information.out.seqsero2_serotype                                                        , remainder: true, by: 0)
  //   .join(contig_information.out.seqsero2_contamination                                                   , remainder: true, by: 0)
  //   .join(contig_information.out.serotypefinder_ogroup                                                    , remainder: true, by: 0)
  //   .join(contig_information.out.serotypefinder_hgroup                                                    , remainder: true, by: 0)
  //   .join(contig_information.out.kraken2_top_hit                                                          , remainder: true, by: 0)
  //   .join(contig_information.out.kraken2_top_perc                                                         , remainder: true, by: 0)
  //   .join(contig_information.out.kraken2_top_reads                                                        , remainder: true, by: 0)
  //   .join(contig_information.out.plasmidfinder_hits                                                       , remainder: true, by: 0)
  //   .join(contig_information.out.quast_gc                                                                 , remainder: true, by: 0)
  //   .join(contig_information.out.quast_contigs                                                            , remainder: true, by: 0)
  //   .join(contig_information.out.quast_nfifty                                                             , remainder: true, by: 0)
  //   .join(contig_information.out.quast_length                                                             , remainder: true, by: 0)
  //   .join(contig_information.out.kleborate_score                                                          , remainder: true, by: 0)
  //   .join(contig_information.out.kleborate_mlst                                                           , remainder: true, by: 0)
  //   .join(contig_information.out.amrfinder_amr_genes                                                      , remainder: true, by: 0)
  //   .join(contig_information.out.amrfinder_vir_genes                                                      , remainder: true, by: 0)
  //   .join(contig_information.out.fastani_ref                                                              , remainder: true, by: 0)
  //   .join(contig_information.out.fastani_ani_score                                                        , remainder: true, by: 0)
  //   .join(contig_information.out.fastani_fragment                                                         , remainder: true, by: 0)
  //   .join(contig_information.out.fastani_total                                                            , remainder: true, by: 0)
  //   .join(contig_information.out.mlst_sttype                                                              , remainder: true, by: 0)

  //   // blobtools
  //   .join(blobtools_species                                                                               , remainder: true, by: 0)
  //   .join(blobtools_perc                                                                                  , remainder: true, by: 0)

  //   // whew!
  //   .set { results }
  // summary(results)

  // contig_information.out.seqsero2_collect
  //   .collectFile(name: "SeqSero_result.tsv",
  //     keepHeader: true,
  //     sort: true,
  //     storeDir: "${params.outdir}/seqsero2")

  // summary.out.summary_files_txt
  //   .collectFile(name: "grandeur_summary.txt",
  //     keepHeader: true,
  //     sort: true,
  //     storeDir: "${params.outdir}/summary")

  // summary.out.summary_files_tsv
  //   .collectFile(name: "grandeur_results.tsv",
  //     keepHeader: true,
  //     sort: true,
  //     storeDir: "${params.outdir}")
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Final Steps

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("MultiQC report can be found at ${params.outdir}/multiqc/multiqc_report.html")
    println("Summary can be found at ${params.outdir}/grandeur_results.tsv")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
