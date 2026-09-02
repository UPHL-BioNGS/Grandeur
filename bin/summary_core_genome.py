import pandas as pd


from bin.summary_warnings import append_warning

# core genome analysis file is also from multiqc
def summarize_core_genome(summary_df, core):
    file = core
    analysis = "core_genome_genes"
    print("Adding core genome percentage from " + file)
    new_df = pd.read_table(file, dtype = str, index_col= False)
    new_df['Sample'] = new_df['Sample'].astype(str)
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_Sample", how = 'left')
    summary_df.drop(analysis + "_Sample", axis=1, inplace=True)
    summary_df['per_core_genome_genes'] = summary_df[analysis + '_core'].astype(float) / (summary_df[analysis + '_soft'].astype(float) + summary_df[analysis + '_core'].astype(float) + summary_df[analysis + '_shell'].astype(float) + summary_df[analysis + '_cloud'].astype(float))
    summary_df['per_core_genome_genes'] = summary_df['per_core_genome_genes'].astype(float) * 100
    summary_df['per_core_genome_genes'] = summary_df['per_core_genome_genes'].round(2)

    return summary_df