import pandas as pd

from bin.summary_warnings import append_warning

# skani
def summarize_skani(summary_df, skani):
    print("Adding results for " + skani)
    analysis = "skani"
    new_df = pd.read_table(skani, dtype=str, index_col=False)
    new_df['sample'] = new_df['sample'].astype(str)
    
    # Clean up names and split Genus_species
    new_df['organism'] = new_df['Ref_file'].str.split('_').str[0:2].str.join('_')
    new_df = new_df.sort_values(['sample', 'ANI'], ascending=[True, False])
    new_df = new_df.add_prefix('skani_')

    top_df = new_df.drop_duplicates(subset=['skani_sample'], keep='first').copy()
    
    # Count how many unique organisms Skani thinks are in this one "isolate"
    counts_df = new_df.groupby('skani_sample').agg(
        skani_predicted_organisms=('skani_organism', 'nunique')
    ).reset_index()
    summary_df = pd.merge(summary_df, top_df, left_on="sample", right_on="skani_sample", how='left').drop('skani_sample', axis=1)
    summary_df = pd.merge(summary_df, counts_df, left_on="sample", right_on="skani_sample", how='left')
    return summary_df