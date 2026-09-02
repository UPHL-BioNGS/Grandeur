import pandas as pd
import numpy as np


from bin.summary_warnings import append_warning

def add_quast_warnings(df):

    if 'quast_n50' in df.columns:
        n50 = pd.to_numeric(df['quast_n50'], errors='coerce')
        append_warning(
            df,
            n50 < 30000,
            "Low N50 (<30000)"
        )

    if 'quast_#_contigs' in df.columns:
        contigs = pd.to_numeric(df['quast_#_contigs'], errors='coerce')
        append_warning(
            df,
            contigs > 500,
            "High contig count (>500)"
        )

    if 'quast_mapped_(%)' in df.columns:
        mapped = pd.to_numeric(df['quast_mapped_(%)'], errors='coerce')
        append_warning(
            df,
            mapped < 90.0,
            "Low mapping rate (<90%)"
        )

    return df

def summarize_quast_from_reads(quast):
    print("Adding results for " + quast)
    file = quast
    analysis = str(file).split("_")[0]
    df = pd.read_table(file, dtype = str, index_col= False)
    df = df.add_prefix(analysis + "_")
    df.columns = [x.lower() for x in df.columns]
    return df

def summarize_quast_contig(quast_contig):
    print("Adding results for " + quast_contig)
    file = quast_contig
    analysis = str(file).split("_")[0]
    df = pd.read_table(file, dtype = str, index_col= False)
    df = df.add_prefix(analysis + "_")
    df.columns = [x.lower() for x in df.columns]
    return df

def summarize_quast(summary_df, q_df, qc_df):
    analysis = "quast"
    new_df = pd.concat([q_df, qc_df])
    new_df[analysis + "_sample"] = new_df[analysis + "_sample"].astype(str)

    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how='left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

    return summary_df