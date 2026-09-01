import pandas as pd

# serotypefinder : splitting O and H groups, getting the top hit for O and H group, combining rows
def summarize_serotypefinder(summary_df, serotypefinder):
    file = serotypefinder
    print("Adding results for " + file)
    analysis = "serotypefinder"
    new_df = pd.read_table(file, dtype=str, index_col=False)
    new_df['sample'] = new_df['sample'].astype(str)
    
    counts_df = new_df.groupby(['sample', 'Database']).size().unstack(fill_value=0).reset_index()
    counts_df = counts_df.rename_axis(None, axis=1) # Clean up the column grouping name
    
    if 'O_type' not in counts_df.columns:
        counts_df['O_type'] = 0
    if 'H_type' not in counts_df.columns:
        counts_df['H_type'] = 0
        
    counts_df = counts_df.rename(columns={
        'O_type': analysis + '_O_count', 
        'H_type': analysis + '_H_count'
    })

    new_df = new_df.sort_values(by='Identity', ascending=False)
    new_df = new_df.drop_duplicates(subset=['sample', 'Database'], keep="first")
    new_df = new_df.add_prefix(analysis + '_')
    
    H_df = new_df[new_df[analysis + '_Database' ] == 'H_type'].copy()
    H_df = H_df.add_suffix('_H')
    
    O_df = new_df[new_df[analysis + '_Database' ] == 'O_type'].copy()
    O_df = O_df.add_suffix('_O')
    
    summary_df = pd.merge(summary_df, O_df, left_on="sample", right_on=analysis + "_sample_O", how='left')
    summary_df.drop(analysis + "_sample_O", axis=1, inplace=True)
    
    summary_df = pd.merge(summary_df, H_df, left_on="sample", right_on=analysis + "_sample_H", how='left')
    summary_df.drop(analysis + "_sample_H", axis=1, inplace=True)
    
    summary_df = pd.merge(summary_df, counts_df[['sample', analysis + '_O_count', analysis + '_H_count']], on="sample", how="left")

    return serotypefinder
