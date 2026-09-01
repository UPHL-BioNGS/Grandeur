import pandas as pd
# (Assuming summary_df and fastqc variables are defined)
def summarize_fastqc(summary_df, fastqc):
    file = fastqc
    print("Adding results for " + file)
    analysis = "fastqc"
    new_df = pd.read_csv(file, dtype=str, index_col=False)
    new_df['sample'] = new_df['sample'].astype(str)
    
    # FastQC output usually has two rows per sample (R1 and R2).
    R1_df = new_df.drop_duplicates(subset='sample', keep="first").add_prefix('R1_')
    R2_df = new_df.drop_duplicates(subset='sample', keep="last").add_prefix('R2_')
    
    new_df = pd.merge(R1_df, R2_df, left_on="R1_sample", right_on="R2_sample", how='left')
    new_df['sample'] = new_df['R1_sample']
    
    # Safely convert to integers
    new_df['R1_Total Sequences'] = new_df['R1_Total Sequences'].astype(int)
    new_df['R2_Total Sequences'] = new_df['R2_Total Sequences'].astype(int)
    new_df['R1_Sequences flagged as poor quality'] = new_df['R1_Sequences flagged as poor quality'].astype(int)
    new_df['R2_Sequences flagged as poor quality'] = new_df['R2_Sequences flagged as poor quality'].astype(int)
    
    # Calculate total and flagged sequences (with the fastqc_ prefix your final script expects)
    new_df['fastqc_total_sequences'] = new_df['R1_Total Sequences'] + new_df['R2_Total Sequences']
    new_df['fastqc_flagged_sequences'] = new_df['R1_Sequences flagged as poor quality'] + new_df['R2_Sequences flagged as poor quality']
    new_df['fastqc_percent_flagged'] = (new_df['fastqc_flagged_sequences'] / new_df['fastqc_total_sequences']) * 100

    # Split ranges like "35-251" by '-' and grab the last item (the max length), then convert to int
    new_df['R1_max_len'] = new_df['R1_Sequence length'].astype(str).apply(lambda x: int(x.split('-')[-1]))
    new_df['R2_max_len'] = new_df['R2_Sequence length'].astype(str).apply(lambda x: int(x.split('-')[-1]))
    
    # Calculate the average read length between R1 and R2
    new_df['fastqc_avg_length'] = (new_df['R1_max_len'] + new_df['R2_max_len']) / 2
    
    # Drop the intermediate R1/R2 calculation columns so they don't clutter
    cols_to_keep = ['sample', 'fastqc_total_sequences', 'fastqc_flagged_sequences', 'fastqc_percent_flagged', 'fastqc_avg_length']
    fastqc_final_df = new_df[cols_to_keep].copy()
    
    # Merge into your main summary dataframe
    summary_df = pd.merge(summary_df, fastqc_final_df, on="sample", how='left')

    return summary_df