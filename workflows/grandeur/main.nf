include { AVERAGE_NUCLEOTIDE_IDENTITY }   from "../../subworkflows/local/average_nucleotide_identity"
include { DE_NOVO_ALIGNMENT }             from "../../subworkflows/local/de_novo_alignment" 
include { PHYLOGENETIC_ANALYSIS }         from "../../subworkflows/local/phylogenetic_analysis"
include { QUALITY_ASSESSMENT }            from "../../subworkflows/local/quality_assessment"
include { REPORT }                        from "../../subworkflows/local/report"
include { SUBTYPING }                     from "../../subworkflows/local/subtyping"
include { TAXONOMIC_PROFILING }           from "../../subworkflows/local/taxonomic_profiling"

workflow GRANDEUR {
    take:
    ch_raw_reads
    ch_fastas
    ch_reference_genomes
    ch_versions
    ch_genome_sizes
    ch_mash_db
    ch_kraken2_db
    ch_checkm2_db
    ch_sylph_db
    dataset_script
    evaluat_script
    jsoncon_script
    multiqc_script
    summary_script
    summfle_script
    version_script

    main:
    ch_for_multiqc   = channel.empty()
    ch_for_summary   = ch_genome_sizes
    ch_for_flag      = channel.empty()
    ch_versions      = channel.empty()
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
        ch_clean_reads   = channel.empty()
        ch_assembled     = channel.empty()
    }

    // getting a summary of everything
    if ( ! params.skip_extras ) {
        // optional subworkflow kraken2 (useful for interspecies contamination)
        TAXONOMIC_PROFILING(
            ch_clean_reads.ifEmpty([]), 
            ch_fastas.ifEmpty([]), 
            ch_assembled.ifEmpty([]), 
            ch_kraken2_db.ifEmpty([]),
            ch_mash_db.ifEmpty([]),
            ch_sylph_db.ifEmpty([])
            )

        ch_for_multiqc = ch_for_multiqc.mix(TAXONOMIC_PROFILING.out.for_multiqc)
        ch_for_summary = ch_for_summary.mix(TAXONOMIC_PROFILING.out.for_summary)
        ch_for_flag    = ch_for_flag.mix(TAXONOMIC_PROFILING.out.for_ref_download)
        ch_versions    = ch_versions.mix(TAXONOMIC_PROFILING.out.versions)

        // determining organisms in sample
        AVERAGE_NUCLEOTIDE_IDENTITY(
            ch_contigs,
            ch_reference_genomes.ifEmpty([]),
            TAXONOMIC_PROFILING.out.for_ref_download.ifEmpty([]),
            dataset_script)

        ch_versions    = ch_versions.mix(AVERAGE_NUCLEOTIDE_IDENTITY.out.versions)
        ch_for_summary = ch_for_summary.mix(AVERAGE_NUCLEOTIDE_IDENTITY.out.for_summary)
        ch_top_hit     = AVERAGE_NUCLEOTIDE_IDENTITY.out.top_hit
        ch_org_contigs = AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_org_contigs

        QUALITY_ASSESSMENT(
            ch_raw_reads.ifEmpty([]),
            ch_clean_reads.ifEmpty([]),
            ch_fastas.ifEmpty([]),
            ch_contigs.ifEmpty([]),
            ch_reads_contigs.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_org_contigs.ifEmpty([]),
            ch_checkm2_db.ifEmpty([]),
            summfle_script)

        ch_for_multiqc = ch_for_multiqc.mix(QUALITY_ASSESSMENT.out.for_multiqc)
        ch_for_summary = ch_for_summary.mix(QUALITY_ASSESSMENT.out.for_summary)
        ch_versions    = ch_versions.mix(QUALITY_ASSESSMENT.out.versions)


        // getting all the other information
        SUBTYPING(
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_myco.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_gas.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_kleb.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_legionella.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_strep.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_salmonella.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_ecoli.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_vibrio.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_gc.ifEmpty([]),
            AVERAGE_NUCLEOTIDE_IDENTITY.out.ch_acinetobacter.ifEmpty([]),
            summfle_script,
            jsoncon_script)

        ch_for_summary = ch_for_summary.mix(SUBTYPING.out.for_summary)
        ch_versions    = ch_versions.mix(SUBTYPING.out.versions)
    } else {
        ch_top_hit = channel.empty()
        ch_org_contigs = ch_contigs.map{it -> tuple(it[0], ["Unknown", "Unknown"], it[1])}
    }

    // optional subworkflow for comparing shared genes
    if ( params.msa ) {
        PHYLOGENETIC_ANALYSIS(
            evaluat_script,
            ch_org_contigs,
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