include { AVERAGE_NUCLEOTIDE_IDENTITY }   from "../subworkflows/local/average_nucleotide_identity"
include { BLOBTOOLS }                     from "../subworkflows/local/blobtools"
include { DE_NOVO_ALIGNMENT }             from "../subworkflows/local/de_novo_alignment" 
include { INFO }                          from "../subworkflows/local/info"
include { KMER_TAXONOMIC_CLASSIFICATION } from "../subworkflows/local/kmer_taxonomic_classification"
include { MIN_HASH }                      from "../subworkflows/local/min_hash"
include { PHYLOGENETIC_ANALYSIS }         from "../subworkflows/local/phylogenetic_analysis"
include { QUALITY_ASSESSMENT }            from "../subworkflows/local/quality_assessment"
include { REPORT }                        from "../subworkflows/local/report"

workflow GRANDEUR {
    take:
    ch_raw_reads
    ch_fastas
    ch_fastani_genomes
    ch_versions
    ch_genome_sizes
    ch_mash_db
    ch_kraken2_db
    ch_blast_db
    dataset_script
    evaluat_script
    jsoncon_script
    multiqc_script
    summary_script
    summfle_script
    version_script

    main:
    ch_for_multiqc   = Channel.empty()
    ch_for_summary   = ch_genome_sizes
    ch_for_flag      = Channel.empty()
    ch_versions      = Channel.empty()
    ch_reads_contigs = ch_fastas.map{it -> tuple(it[0], it[1], null)}


    if ( params.sample_sheet || params.reads || params.sra_accessions ) {
        DE_NOVO_ALIGNMENT(ch_raw_reads)

        ch_assembled     = DE_NOVO_ALIGNMENT.out.contigs
        ch_contigs       = ch_fastas.mix(DE_NOVO_ALIGNMENT.out.contigs)
        ch_reads_contigs = ch_reads_contigs.mix(DE_NOVO_ALIGNMENT.out.reads_contigs)
        ch_clean_reads   = DE_NOVO_ALIGNMENT.out.clean_reads
        ch_for_multiqc   = ch_for_multiqc.mix(DE_NOVO_ALIGNMENT.out.for_multiqc)
        ch_versions      = ch_versions.mix(DE_NOVO_ALIGNMENT.out.versions)

    } else {
        ch_contigs       = ch_fastas
        ch_clean_reads   = Channel.empty()
        ch_assembled     = Channel.empty()
    }

    // getting a summary of everything
    if ( ! params.skip_extras ) {
        QUALITY_ASSESSMENT(
            ch_raw_reads,
            ch_contigs,
            ch_reads_contigs,
            summfle_script)

        ch_for_multiqc = ch_for_multiqc.mix(QUALITY_ASSESSMENT.out.for_multiqc)
        ch_for_summary = ch_for_summary.mix(QUALITY_ASSESSMENT.out.for_summary)
        ch_versions    = ch_versions.mix(QUALITY_ASSESSMENT.out.versions)


        // optional subworkflow blobtools (useful for interspecies contamination)
        if ( params.blast_db && ( params.sample_sheet || params.reads || params.sra_accessions )) {
            BLOBTOOLS(QUALITY_ASSESSMENT.out.bams, ch_blast_db )

            ch_for_summary = ch_for_summary.mix(BLOBTOOLS.out.for_summary)
            ch_for_flag    = ch_for_flag.mix(BLOBTOOLS.out.for_flag)
            ch_versions = ch_versions.mix(BLOBTOOLS.out.versions)
        }

        // optional subworkflow kraken2 (useful for interspecies contamination)
        if ( params.kraken2_db && ( params.sample_sheet || params.reads || params.sra_accessions )) {
            KMER_TAXONOMIC_CLASSIFICATION(ch_clean_reads, ch_kraken2_db )

            ch_for_multiqc = ch_for_multiqc.mix(KMER_TAXONOMIC_CLASSIFICATION.out.for_multiqc)
            ch_for_summary = ch_for_summary.mix(KMER_TAXONOMIC_CLASSIFICATION.out.for_summary)
            ch_for_flag    = ch_for_flag.mix(KMER_TAXONOMIC_CLASSIFICATION.out.for_flag)
            ch_versions    = ch_versions.mix(KMER_TAXONOMIC_CLASSIFICATION.out.versions)
        } 

        // subworkflow mash for species determination
        MIN_HASH(ch_clean_reads, ch_fastas, ch_mash_db)
        ch_versions = ch_versions.mix(MIN_HASH.out.versions)
        ch_for_summary = ch_for_summary.mix(MIN_HASH.out.for_summary)

        // determining organisms in sample
        AVERAGE_NUCLEOTIDE_IDENTITY(
            ch_for_summary.collect(),
            ch_contigs,
            ch_fastani_genomes.ifEmpty([]),
            dataset_script)

        ch_versions = ch_versions.mix(AVERAGE_NUCLEOTIDE_IDENTITY.out.versions)
        ch_for_flag = ch_for_flag.mix(AVERAGE_NUCLEOTIDE_IDENTITY.out.for_flag).mix(MIN_HASH.out.for_flag)
        ch_top_hit  = AVERAGE_NUCLEOTIDE_IDENTITY.out.top_hit
        ch_for_summary = ch_for_summary.mix(AVERAGE_NUCLEOTIDE_IDENTITY.out.for_summary)


        // getting all the other information
        INFO(
            ch_contigs, 
            ch_for_flag, 
            summfle_script,
            jsoncon_script)

        ch_for_summary = ch_for_summary.mix(INFO.out.for_summary)
        ch_versions    = ch_versions.mix(INFO.out.versions)
    } else {
        ch_top_hit = Channel.empty()
    }

    // optional subworkflow for comparing shared genes
    if ( params.msa ) {
        PHYLOGENETIC_ANALYSIS(
            evaluat_script,
            ch_contigs.ifEmpty([]),
            ch_top_hit.ifEmpty([]))
            
        ch_for_multiqc = ch_for_multiqc.mix(PHYLOGENETIC_ANALYSIS.out.for_multiqc)
        ch_versions    = ch_versions.mix(PHYLOGENETIC_ANALYSIS.out.versions)
    }

    // getting a summary of everything
    if ( ! params.skip_extras ) {
        REPORT(
            ch_raw_reads, 
            ch_fastas, 
            ch_for_multiqc.collect(), 
            ch_for_summary.concat(summary_script).collect(),
            ch_versions.collect(),
            multiqc_script,
            version_script
        )
    }
}