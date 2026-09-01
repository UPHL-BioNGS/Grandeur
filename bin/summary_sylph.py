import pandas as pd

# sylph
def summarize_sylph(summary_df, sylph):
    print("Adding results for " + sylph)
    new_df = pd.read_table(sylph, dtype=str, index_col=False)
    new_df['sample'] = new_df['sample'].astype(str)
    
    # Convert abundance to numeric so we can safely sort by it
    new_df['Taxonomic_abundance'] = pd.to_numeric(new_df['Taxonomic_abundance'], errors='coerce')
    
    # Sort to ensure the top hit (highest abundance) is first for every sample
    new_df = new_df.sort_values(['sample', 'Taxonomic_abundance'], ascending=[True, False])
    
    # --- 1. Get the Top Hit ---
    # Keep only the first row per sample to maintain the main dataframe columns
    top_df = new_df.drop_duplicates(subset=['sample'], keep='first').copy()
    top_df = top_df.add_prefix('sylph_')
    
    # --- 2. Count Unique Organisms ---
    # Group by sample and count the unique Genome_files
    counts_df = new_df.groupby('sample').agg(
        sylph_predicted_organisms=('Genome_file', 'nunique')
    ).reset_index()
    
    # --- 3. Merge ---
    # Merge the top hit details
    summary_df = pd.merge(summary_df, top_df, left_on="sample", right_on="sylph_sample", how='left').drop('sylph_sample', axis=1)
    
    # Merge the organism count
    summary_df = pd.merge(summary_df, counts_df, on="sample", how='left')

    return summary_df