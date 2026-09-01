import pandas as pd

# kraken2 : merging relevant rows into one
def summarize_kraken2(summary_df, kraken2):
    print("Adding results for " + kraken2)
    analysis = "kraken2"
    new_df = pd.read_csv(kraken2, dtype=str, index_col=False)
    new_df['Sample'] = new_df['Sample'].astype(str)
    
    # Sort by abundance to find the top hit
    new_df = new_df.sort_values(['Sample', 'Percentage of fragments'], ascending=False)
    
    # Save the top organism name for the "Predicted Organism" logic later
    tmp_df = new_df.drop_duplicates(subset=['Sample'], keep="first")[['Sample', 'Scientific name']]
    tmp_df.columns = ['Sample', 'kraken2_top_organism']

    # Create the string format for the report (e.g., "E.coli (98%)")
    new_df['kraken2_formatted'] = new_df['Scientific name'] + " (" + new_df['Percentage of fragments'] + "%)"
    
    # Group by Sample and calculate both the report list and the unique organism count
    grouped = new_df.groupby('Sample').agg(
        kraken2_report=('kraken2_formatted', lambda x: ', '.join(x.dropna().astype(str))),
        kraken2_predicted_organisms=('Scientific name', 'nunique')
    ).reset_index()
    
    # Merge both dataframes into the main summary_df
    summary_df = pd.merge(summary_df, tmp_df, left_on="sample", right_on="Sample", how='left').drop('Sample', axis=1)
    summary_df = pd.merge(summary_df, grouped, left_on="sample", right_on="Sample", how='left').drop('Sample', axis=1)

    return summary_df