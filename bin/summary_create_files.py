
def create_final_summary(df):
# reducing to the top 1 or 2 results for each analysis
    final_columns = [
        # general information
        'coverage',
        'per_core_genome_genes',
        'fastqc_total_sequences',
        'fastqc_flagged_sequences',
        'fastqc_avg_length',
        'fastp_pct_passed_reads',
        'quast_#_contigs',
        'quast_mapped_(%)',
        'quast_gc_(%)',
        'checkm2_completeness',
        'checkm2_contamination',
        'warnings',
        'amrfinder_genes_(per_cov/per_ident)',

        # species
        'predicted_organism',
        'mlst_matching_pubmlst_scheme',
        'mlst_st',
        'skani_ANI',
        'skani_organism',
        'spestimator_num_refs',
        'datasets_num_genomes',
        'sylph_Adjusted_ANI',
        'sylph_Taxonomic_abundance',
        'mash_dist_organism',
        'mash_dist_mash-distance',
        'mash_dist_p-value',
        'mash_screen_organism',
        'mash_screen_identity',
        'mash_screen_shared-hashes',
        'mash_screen_p-value',
        'plasmidfinder_plasmid_(identity)',

        # contamination
        'kraken2_report',

        # species specific information
        'seqsero2_predicted_antigenic_profile',
        'seqsero2_predicted_serotype',
        'emmtyper_predicted_emm-type',
        'kleborate_virulence_score',
        'kleborate_resistance_score',
        'kaptive_acinetobacter_baumannii_oc_locus_primary_reference_best_match_locus',
        'kaptive_acinetobacter_baumannii_k_locus_primary_reference_best_match_locus',
        'kaptive_klebsiella_k_locus_primary_reference_best_match_locus',
        'kaptive_klebsiella_o_locus_primary_reference_best_match_locus',
        'kaptive_vibriopara_kaptivedb_k_best_match_locus',
        'kaptive_vibriopara_kaptivedb_o_best_match_locus',
        'elgato_st',
        'meningotype_serogroup',
        'ngmaster_ng-mast/ng-star',
        'mykrobe_phylo_group',
        'mykrobe_species',
        'mykrobe_lineage',
        'drprg_susceptibility',
        'pbptyper_pbptype',
        'serotypefinder_Serotype_O',
        'serotypefinder_Serotype_H',
        'shigapass_ipah',
        'shigapass_predicted_serotype',

        # phylogenetic analysis results
        'snpdists_core_gene_alignment.txt_average',
        'snpdists_core_gene_alignment.txt_min',
        'snpdists_core_gene_alignment.txt_max',
        'iqtree_core_gene_alignment.treefile.nwk_average',
        'iqtree_core_gene_alignment.treefile.nwk_min',
        'iqtree_core_gene_alignment.treefile.nwk_max',
        'iqtree_core_gene_alignment_meanbrlen',
        'iqtree_core_gene_alignment_sumbrlen',
        'snpdists_ska_alignment.txt_average',
        'snpdists_ska_alignment.txt_min',
        'snpdists_ska_alignment.txt_max',
        'iqtree_ska_alignment.treefile.nwk_average',
        'iqtree_ska_alignment.treefile.nwk_min',
        'iqtree_ska_alignment.treefile.nwk_max',
        'iqtree_ska_alignment_meanbrlen',
        'iqtree_ska_alignment_sumbrlen',
        'mashtree.nwk_average',
        'mashtree.nwk_min',
        'mashtree.nwk_max',
        'mashtree_meanbrlen',
        'mashtree_sumbrlen'

        ]

    set_columns = []
    for new_column in final_columns :
        if new_column in summary_df.columns :
            set_columns.append(new_column)

    summary_df.to_csv(final + '.tsv', columns = ['sample','file','version'] + set_columns, index=False, sep="\t")
    summary_df.to_csv(final + '.txt', columns = ['sample','file','version'] + set_columns, index=False, sep=";")

def create_extended_summary(summary_df):

    ##########################################
    # creating files                         #
    ##########################################

    print("Creating final files")

    summary_df = summary_df.sort_values(by='sample')
    summary_df = summary_df.fillna("")

    summary_df.columns = summary_df.columns.str.replace(' ', '_')

    summary_df.to_csv(extended + '.tsv', index=False, sep="\t")
    summary_df.to_csv(extended + '.txt', index=False, sep=";")

    transposed_df = summary_df.set_index('sample').T
    transposed_df = transposed_df.reset_index().rename(columns={'index': 'metric'})
    transposed_df.to_csv(extended + '_transposed.tsv', index=False, sep="\t")
    transposed_df.to_csv(extended + '_transposed.txt', index=False, sep=";")


    