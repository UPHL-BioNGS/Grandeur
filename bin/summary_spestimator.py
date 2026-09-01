import pandas as pd

# spestimator : counting unique reference hits per sample
def summarize_spestimator(summary_df, spestimator):
    print("Adding results for " + spestimator)
    analysis = "spestimator"
    
    # Read the TSV file
    new_df = pd.read_table(spestimator, sep=",", dtype=str)
    
    # 1. Clean the 'input file' column to match your 'sample' IDs
    # This removes '_contigs.fa' and other extensions
    new_df['sample'] = (
        new_df['input file']
        .str.replace('_contigs.fa', '', regex=False)
        .str.replace(r'\.(fasta|fna|fa)$', '', regex=True)
    )
    
    new_df['sample'] = new_df['sample'].astype(str)
    
    # 2. Count unique organisms identified for each sample
    # Using 'organism' as the unique identifier here
    sp_counts = new_df.groupby('sample')['organism'].nunique().reset_index()
    sp_counts.columns = ['sample', 'spestimator_num_refs']
    
    # 3. Merge into the main summary_df
    summary_df = pd.merge(summary_df, sp_counts, on='sample', how='left')

    return summary_df