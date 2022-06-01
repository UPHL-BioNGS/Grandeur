#!/usr/bin/env nextflow

println("Currently using the Grandeur workflow for use with microbial sequencing. The view is great from 8299 feet (2530 meters) above sea level.\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v2.0.20220610")
println("")

nextflow.enable.dsl = 2

params.outdir                     = workflow.launchDir + '/grandeur'
params.maxcpus                    = 12
params.medcpus                    = 4

// core workflow of fastq to contig
params.fastp_options              = "--detect_adapter_for_pe"
params.bbduk_options              = "k=31 hdist=1"
params.spades_options             = '--isolate'

// fastq information
// 'shigatyper',
params.fastq_processes            = ['fastp', 'bbduk', 'spades', 'fastqc', 'cg_pipeline',  'mash', 'kraken2', 'summary', 'seqsero2', 'multiqc']
params.fastqc_options             = ''
params.cg_pipeline_options        = '--qual_offset 33 --minLength 1'
params.shigatyper_options         = ''
params.kraken2_db                 = false
params.kraken2_options            = ''
// WARNING : DO NOT CHANGE params.mash_reference UNLESS YOU ARE USING A DIFFERENT MASH CONTAINER
params.mash_reference             = '/db/RefSeqSketchesDefaults.msh'
params.mash_options               = '-v 0 -d 0.5'
// duplicates as seqsero2 and serotypefinder defaults are for contigs
// params.seqsero2_options        = '-t 2 -m a -b mem'
// params.serotypefinder_options  = ''

// contig information
params.contig_processes           = ['amrfinderplus', 'kleborate', 'fastani', 'mlst', 'quast', 'serotypefinder', 'blobtools', 'summary', 'multiqc']
params.amrfinderplus_options      = ''
params.fastani_options            = ''
params.kleborate_options          = '-all'
params.mlst_options               = ''
params.quast_options              = ''
params.serotypefinder_options     = ''
params.seqsero2_options           = '-m a -b mem'
// // duplicate as kraken2 with fastq reads is default
// // params.kraken2_options         = ''

// for blobtools
params.blast_db                   = false
params.local_db_type              = 'nt'
params.blobtools_create_options   = ''
params.blobtools_view_options     = ''
params.blobtools_plot_options     = '--format png -r species'
params.bwa_options                = ''
params.samtools_sort_options      = ''

// summary
params.multiqc_options            = ''

// phylogenetic analysis
params.phylogenetic_processes     = []
params.iqtree2_options            = '-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000'
params.outgroup                   = ''
params.prokka_options             = '--mincontiglen 500 --compliant --locustag locus_tag --centre STAPHB'
params.roary_options              = ''
params.snp_dists_options          = '-c'

include { de_novo_alignment }     from './subworkflows/de_novo_alignment.nf'     addParams( outdir:                     params.outdir,
                                                                                            fastq_processes:            params.fastq_processes,
                                                                                            spades_options:             params.spades_options,
                                                                                            bbduk_optons:               params.bbduk_options,
                                                                                            fastp_options:              params.fastp_options)
include { fastq_information }     from './subworkflows/fastq_information.nf'     addParams( outdir:                     params.outdir,
                                                                                            fastq_processes:            params.fastq_processes,
                                                                                            fastqc_options:             params.fastqc_options,
                                                                                            cg_pipeline_options:        params.cg_pipeline_options,
                                                                                            shigatyper_options:         params.shigatyper_options,
                                                                                            kraken2_options:            params.kraken2_options,
                                                                                            mash_reference:             params.mash_reference,
                                                                                            mash_options:               params.mash_options,
                                                                                            serotypefinder_options:     params.serotypefinder_options,
                                                                                            seqsero2_options:           params.seqsero2_options)
include { contig_information }    from './subworkflows/contig_information.nf'    addParams( outdir:                     params.outdir,
                                                                                            contig_processes:           params.contig_processes,
                                                                                            amrfinderplus_options:      params.amrfinderplus_options,
                                                                                            fastani_options:            params.fastani_options,
                                                                                            kleborate_options:          params.kleborate_options,
                                                                                            mlst_options:               params.mlst_options,
                                                                                            quast_options:              params.quast_options,
                                                                                            serotypefinder_options:     params.serotypefinder_options,
                                                                                            seqsero2_options:           params.seqsero2_options,
                                                                                            kraken2_options:            params.kraken2_options)
include { blobtools }             from './subworkflows/blobtools.nf'             addParams( local_db_type:              params.local_db_type,
                                                                                            blobtools_create_options:   params.blobtools_create_options,
                                                                                            blobtools_view_options:     params.blobtools_view_options,
                                                                                            blobtools_plot_options:     params.blobtools_plot_options,
                                                                                            bwa_options:                params.bwa_options,
                                                                                            samtools_sort_options:      params.samtools_sort_options)
include { phylogenetic_analysis } from './subworkflows/phylogenetic_analysis.nf' addParams( outdir:                     params.outdir,
                                                                                            phylogenetic_processes:     params.phylogenetic_processes,
                                                                                            iqtree2_options:            params.iqtree2_options,
                                                                                            outgroup:                   params.outgroup,
                                                                                            prokka_options:             params.prokka_options,
                                                                                            roary_options:              params.roary_options,
                                                                                            snp_dists_options:          params.snp_dists_options)
include { mash_dist as mash }     from './modules/mash'                          addParams( fastq_processes:            params.fastq_processes,
                                                                                            mash_reference:             params.mash_reference,
                                                                                            mash_options:               params.mash_options)
include { summary }               from './modules/summary'                       addParams( fastq_processes:            params.fastq_processes,
                                                                                            contig_processes:           params.contig_processes,
                                                                                            phylogenetic_processes:     params.phylogenetic_processes)
include { multiqc }               from './modules/multiqc'                       addParams( multiqc_options:            params.multiqc_options
                                                                                            fastq_processes:            params.fastq_processes,
                                                                                            contig_processes:           params.contig_processes,
                                                                                            phylogenetic_processes:     params.phylogenetic_processes)

// TODO : something for plasmids
// TODO : frp_plasmid
// TODO : pointfinder
// TODO : mubsuite
// TODO : sistr
// TODO : plasmidseeker
// TODO : socru?

// Getting the fastq files
params.reads = workflow.launchDir + '/reads'
Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                  "${params.reads}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .view { "Paired-end fastq files found : ${it[0]}" }
  .set { reads }

// input_reads
//
//   .view { "Paired-end fastq files found : ${it[0]}" }
//   .ifEmpty{
//     println("No fastq or fastq.gz files were found at ${params.reads}")
//     println("Set 'params.reads' to directory with paired-end reads" )
//   }
//   .set { reads }

// Getting contig or fasta files
params.fastas = workflow.launchDir + '/fastas'
Channel
  .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
  .map { file -> tuple(file.baseName, file) }
  .set { fastas }

// input_fastas
//   .map { file -> tuple(file.baseName, file) }
//   .view { "Fasta file found : ${it[0]}" }
//   .into { fastas_check ; fastas ; fastas_mash ; fastas_quast ; fastas_prokka ; fastas_seqsero2 ; fastas_amrfinder ; fastas_serotypefinder ; fastas_kleborate ; fastas_mlst ; fastas_kraken2 ; fastas_fastani }
//

// Getting fasta files that have been annotated with prokka
params.gff = workflow.launchDir + '/gff'
Channel.fromPath("${params.gff}/*.gff", type: 'file')
  .view { "gff file : $it" }
  .set { gffs }
//
// reads_check
//   .concat(fastas_check)
//   .concat(gffs_check)
//   .ifEmpty{
//     println("FATAL : No fastq or fastq.gz files were found at ${params.reads}")
//     println("Set 'params.reads' to directory with paired-end reads" )
//     println("FATAL : No fastas were found at ${params.fastas}")
//     println("Set 'params.fastas' to directory with fastas" )
//     println("FATAL : No gff files were found at ${params.gff}")
//     println("Set 'params.gff' to directory with gff files" )
//     exit 1
//   }
//
// Getting the file with genome sizes of common organisms for cg-pipeline. The End User can use their own file and set with a param
params.genome_sizes = workflow.projectDir + "/configs/genome_sizes.json"
genome_sizes = Channel.fromPath(params.genome_sizes, type:'file')
//
// Getting the reference genomes for fastANI
params.fastani_refs = workflow.projectDir + "/configs/fastani_ref.tar.gz"
fastani_genomes = Channel.fromPath("${params.fastani_refs}", type:'file')

// Getting the database for blobtools
local_blastdb = params.blast_db
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
local_kraken2 = params.kraken2_db
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

workflow {
  de_novo_alignment(reads)

  fastq_information(reads, de_novo_alignment.out.clean_reads, genome_sizes, local_kraken2)

  if (params.blast_db) {
    blobtools(de_novo_alignment.out.clean_reads, de_novo_alignment.out.contigs, local_blastdb)
    blobtools_species = blobtools.out.species
    blobtools_perc    = blobtools.out.perc
  } else {
    blobtools_species = Channel.empty()
    blobtools_perc    = Channel.empty()
  }

  contigs             = de_novo_alignment.out.contigs.mix(fastas)
  mash(fastas)
  mash_species        = fastq_information.out.mash_species.mix(mash.out.species)
  mash_genus          = fastq_information.out.mash_genus.mix(mash.out.genus)
  salmonella_flag     = fastq_information.out.salmonella_flag.mix(mash.out.salmonella_flag)
  ecoli_flag          = fastq_information.out.ecoli_flag.mix(mash.out.ecoli_flag)
  klebsiella_flag     = fastq_information.out.klebsiella_flag.mix(mash.out.klebsiella_flag)

  contig_information(contigs, mash_species, mash_genus, salmonella_flag, ecoli_flag, klebsiella_flag, fastani_genomes, local_kraken2)

  phylogenetic_analysis(contigs, mash_species, mash_genus, gffs )

  multiqc(de_novo_alignment.out.fastp_multiqc.collect().ifEmpty([]),
          de_novo_alignment.out.bbduk_multiqc.collect().ifEmpty([]),
          contig_information.out.kraken2_multiqc.collect().ifEmpty([]),
          contig_information.out.quast_multiqc.collect().ifEmpty([]),
          fastq_information.out.fastqc_multiqc.collect().ifEmpty([]),
          fastq_information.out.kraken2_multiqc.collect().ifEmpty([]),
          phylogenetic_analysis.out.prokka_multiqc.collect().ifEmpty([]))

  reads
    .mix(fastas)
    // de_novo_alignment
    .join(de_novo_alignment.out.phix_reads                                                                , remainder: true, by: 0)
    .join(de_novo_alignment.out.fastp_results                                                             , remainder: true, by: 0)

    // fastq_information
    .join(fastq_information.out.fastqc_1_results                                                          , remainder: true, by: 0)
    .join(fastq_information.out.fastqc_2_results                                                          , remainder: true, by: 0)
    .join(fastq_information.out.cg_pipeline_read_length                                                   , remainder: true, by: 0)
    .join(fastq_information.out.cg_pipeline_quality                                                       , remainder: true, by: 0)
    .join(fastq_information.out.cg_pipeline_coverage                                                      , remainder: true, by: 0)
    .join(fastq_information.out.cg_pipeline_ref_gen_len                                                   , remainder: true, by: 0)
    .join(fastq_information.out.shigatyper_predictions                                                    , remainder: true, by: 0)
    .join(fastq_information.out.shigatyper_cada                                                           , remainder: true, by: 0)
    .join(fastq_information.out.mash_genome_size                                                          , remainder: true, by: 0)
    .join(fastq_information.out.mash_coverage                                                             , remainder: true, by: 0)

    // fastq_information or contig_information
    .join(mash_genus                                                                                      , remainder: true, by: 0)
    .join(mash_species                                                                                    , remainder: true, by: 0)
    .join(fastq_information.out.mash_full.mix(mash.out.full)                                              , remainder: true, by: 0)
    .join(fastq_information.out.mash_pvalue.mix(mash.out.pvalue)                                          , remainder: true, by: 0)
    .join(fastq_information.out.mash_distance.mix(mash.out.distance)                                      , remainder: true, by: 0)
    .join(fastq_information.out.seqsero2_profile.mix(contig_information.out.seqsero2_profile)             , remainder: true, by: 0)
    .join(fastq_information.out.seqsero2_serotype.mix(contig_information.out.seqsero2_serotype)           , remainder: true, by: 0)
    .join(fastq_information.out.seqsero2_contamination.mix(contig_information.out.seqsero2_contamination) , remainder: true, by: 0)
    .join(fastq_information.out.serotypefinder_ogroup.mix(contig_information.out.serotypefinder_ogroup)   , remainder: true, by: 0)
    .join(fastq_information.out.serotypefinder_hgroup.mix(contig_information.out.serotypefinder_hgroup)   , remainder: true, by: 0)
    .join(fastq_information.out.kraken2_top_hit.mix(contig_information.out.kraken2_top_hit)               , remainder: true, by: 0)
    .join(fastq_information.out.kraken2_top_perc.mix(contig_information.out.kraken2_top_perc)             , remainder: true, by: 0)
    .join(fastq_information.out.kraken2_top_reads.mix(contig_information.out.kraken2_top_reads)           , remainder: true, by: 0)

    // contig_information
    .join(contig_information.out.quast_gc                                                                 , remainder: true, by: 0)
    .join(contig_information.out.quast_contigs                                                            , remainder: true, by: 0)
    .join(contig_information.out.quast_nfifty                                                             , remainder: true, by: 0)
    .join(contig_information.out.quast_length                                                             , remainder: true, by: 0)
    .join(contig_information.out.kleborate_score                                                          , remainder: true, by: 0)
    .join(contig_information.out.kleborate_mlst                                                           , remainder: true, by: 0)
    .join(contig_information.out.amrfinder_amr_genes                                                      , remainder: true, by: 0)
    .join(contig_information.out.amrfinder_vir_genes                                                      , remainder: true, by: 0)
    .join(contig_information.out.fastani_ref                                                              , remainder: true, by: 0)
    .join(contig_information.out.fastani_ani_score                                                        , remainder: true, by: 0)
    .join(contig_information.out.fastani_fragment                                                         , remainder: true, by: 0)
    .join(contig_information.out.fastani_total                                                            , remainder: true, by: 0)
    .join(contig_information.out.mlst_sttype                                                              , remainder: true, by: 0)

    // blobtools
    .join(blobtools_species                                                                               , remainder: true, by: 0)
    .join(blobtools_perc                                                                                  , remainder: true, by: 0)

    // whew!
    .set { results }
  summary(results)

  fastq_information.out.seqsero2_collect
    .concat(contig_information.out.seqsero2_collect)
    .collectFile(name: "SeqSero_result.tsv",
      keepHeader: true,
      sort: true,
      storeDir: "${params.outdir}/seqsero2")

  summary.out.summary_files_txt
    .collectFile(name: "grandeur_summary.txt",
      keepHeader: true,
      sort: true,
      storeDir: "${params.outdir}/summary")

  summary.out.summary_files_tsv
    .collectFile(name: "grandeur_results.tsv",
      keepHeader: true,
      sort: true,
      storeDir: "${params.outdir}")
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("MultiQC report can be found at ${params.outdir}/multiqc/multiqc_report.html")
    println("Summary can be found at ${params.outdir}/grandeur_results.tsv")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
