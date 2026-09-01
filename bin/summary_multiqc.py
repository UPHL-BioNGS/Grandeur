import pandas as pd

if exists(multiqc_stats) : 
    file = multiqc_stats
    print("Adding analysis parsed via multiqc in " + file)
    new_df = pd.read_table(file, dtype = str, index_col= False)
    fastp_columns = [col for col in new_df.columns if col.startswith('fastp')]

    if fastp_columns:
        tmp_df = new_df[["Sample"] + fastp_columns].copy()
        tmp_df["Sample"] = tmp_df["Sample"].astype(str)
        tmp_df["possible_fastp_name"] = tmp_df['Sample'].str.split(" ").str[0].str.split(".").str[0].str.split("_").str[0]
        if 'fastp-pct_surviving' in tmp_df.columns:
            tmp_df = tmp_df.dropna(subset=['fastp-pct_surviving'])
            tmp_df["fastp_pct_passed_reads"] = tmp_df["fastp-pct_surviving"].astype(float).round(2)
            tmp_df.drop("fastp-pct_surviving", axis=1, inplace=True)
        
        summary_df["possible_fastp_name"] = summary_df['file'].str.split(" ").str[0].str.split(".").str[0].str.split("_").str[0]
        summary_df = pd.merge(summary_df, tmp_df, on="possible_fastp_name", how = 'left')
        summary_df.drop("Sample", axis=1, inplace=True)
        summary_df.drop("possible_fastp_name", axis=1, inplace=True)