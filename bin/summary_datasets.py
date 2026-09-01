import pandas as pd

# datasets : adding a count of reference genomes available
def summary_datasets(summary_df, datasets):
    print("Adding reference genome counts from " + datasets)
    # Read the datasets file
    datasets_df = pd.read_csv(datasets, dtype=str)
    
    # Calculate the total number of genomes in the database
    num_genomes = len(datasets_df)
    
    # Assign this constant value to a new column in the main summary
    summary_df['datasets_num_genomes'] = num_genomes

    return summary_df