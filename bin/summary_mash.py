import pandas as pd

# mash dist
def summarize_mashdist(summary_df, mash_dist):
    file = mash_dist
    print("Adding results for " + file)
    analysis = "mash_dist"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    new_df['sample'] = new_df['sample'].astype(str)
    # header : sample,reference,query,mash-distance,P-value,matching-hashes,organism
    new_df = new_df.sort_values(by = ['sample', 'P-value', 'mash-distance'], ascending = [True, True, True])

    counts = new_df.groupby('sample')['organism'].nunique().reset_index()
    counts.columns = ['sample', 'predicted_organisms']

    new_df = pd.merge(new_df, counts, on='sample', how='left')
    new_df = new_df.drop_duplicates(subset=['sample'], keep='first').reset_index(drop=True)
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

    return summary_df

# mash screen
def summarize_mashscreen(summary_df, mash_screen):
    file = mash_screen
    print("Adding results for " + file)
    analysis = "mash_screen"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    new_df['sample'] = new_df['sample'].astype(str)
    # header : sample,identity,shared-hashes,median-multiplicity,p-value,query-ID,organism
    new_df = new_df.sort_values(by = ['sample', 'p-value', 'identity'], ascending = [True, True, False])

    counts = new_df.groupby('sample')['organism'].nunique().reset_index()
    counts.columns = ['sample', 'predicted_organisms']

    new_df = new_df.drop_duplicates(subset=['sample'], keep='first')
    new_df = pd.merge(new_df, counts, on='sample', how='left')
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

    return summary_df