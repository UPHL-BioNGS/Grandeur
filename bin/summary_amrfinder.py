import pandas as pd

# Create the detailed header "Type_Family_Gene"
def element_family(row):
    if 'family' in str(row['Element name']).lower():
        family = str(row['Element name']).lower().split(' family')[0] + " family"
        return f"{row['Type']}_{family}_{row['Element symbol']}"
    return f"{row['Type']}_{row['Element symbol']}"

def summarize_amrfinder(df, file):
    # takes the results from amrfinderplus and summarizes them into a single row per sample in the summary_df
    # the row that goes in the final summary file has one column for the list of genes
    # the standalone matrix file will have one column per gene and gene family with the value being "Detected" | "Not Detected" | "Paritial" | "Point"
    print("Adding results for " + file)
    analysis = "amrfinder"
    new_df = pd.read_table(file, dtype=str, index_col=False)

    # creating a new column with a list of genes and their coverage/identity values
    summary_prep = new_df.copy()
    summary_prep = summary_prep.sort_values('Element symbol')
    summary_prep['genes_formatted'] = summary_prep['Element symbol'] + ' (' + \
        summary_prep['% Coverage of reference'] + '/' + \
        summary_prep['% Identity to reference'] + ')'
        
    # Group by Sample (Name) and make a list
    summary_genes = summary_prep.groupby('Name', as_index=False).agg(
        {'genes_formatted': lambda x: ', '.join(x.dropna().astype(str))}
    )
    summary_genes.columns = ['sample', 'amrfinder_genes_(per_cov/per_ident)']
        
    # Merge into summary_df
    summary_df = pd.merge(df, summary_genes, on="sample", how='left')


    # TODO : create a massive matrix with one column per gene
    # amr_df = new_df.copy()    
    # amr_df['gene_header'] = amr_df.apply(element_family, axis=1)
    # TODO : finish creating AMR family csv file
    # Pivot to create the wide matrix
    # amr_matrix = amr_df.pivot(index='Name', columns='gene_header', values='matrix_val').fillna('-')
        
    # Save the matrix to a separate file
    # amr_matrix.to_csv('amr_virulence_matrix.tsv', sep='\t')
    # print("Standalone AMR matrix saved to amr_virulence_matrix.tsv")

    return summary_df