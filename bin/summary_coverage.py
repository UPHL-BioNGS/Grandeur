import pandas as pd
import numpy as np

##########################################
# size and coverage estimates            #
##########################################
from bin.summary_warnings import append_warning

def add_warnings(summary_df):

    ##########################################
    # warnings and flags                     #
    ##########################################

    print("Adding warnings for QC")

    # coverage flags
    mask_under_30 = summary_df['coverage'] < 30
    mask_under_40 = (summary_df['coverage'] >= 30) & (summary_df['coverage'] < 40)

    # Define the warning messages
    warn_30 = "Low coverage (<30x)"
    warn_40 = "Low coverage (<40x)"

    # Apply the warnings using numpy to conditionally add a semicolon separator if needed
    summary_df.loc[mask_under_30, 'warnings'] += np.where(
        summary_df.loc[mask_under_30, 'warnings'] == '', 
        warn_30, 
        ', ' + warn_30
    )

    summary_df.loc[mask_under_40, 'warnings'] += np.where(
        summary_df.loc[mask_under_40, 'warnings'] == '', 
        warn_40, 
        ', ' + warn_40
    )


def summarize_coverage(summary_df, genome_sizes):

    # 1. Initialize final_coverage
    summary_df['final_coverage'] = np.nan

    # 2. Identify paired-end / fastq rows
    if 'file_2' in summary_df.columns:
        is_paired_end = summary_df['file_2'].notna() & (summary_df['file_2'] != '')
    else:
        # Failsafe: if the column is entirely missing, treat as False
        is_paired_end = pd.Series(False, index=summary_df.index) 

    # 3. Calculate Coverage
    if "fastqc_total_sequences" in summary_df.columns and 'fastqc_avg_length' in summary_df.columns:
        print("Estimating coverage")
        
        total_seqs = pd.to_numeric(summary_df['fastqc_total_sequences'], errors='coerce')
        avg_len = pd.to_numeric(summary_df['fastqc_avg_length'], errors='coerce')
        summary_df['total_bases'] = total_seqs * avg_len
        
        if 'genome_sizes' in locals() and exists(genome_sizes):
            with open(genome_sizes, 'r') as f:
                data = json.load(f)
                size_dict = data.get("genome_sizes", {})
            summary_df['expected_size'] = summary_df['predicted_organism'].map(size_dict)
        else:
            summary_df['expected_size'] = np.nan
            
        expected_size_num = pd.to_numeric(summary_df['expected_size'], errors='coerce')
        summary_df['cov_expected'] = summary_df['total_bases'] / expected_size_num
        
        quast_len = pd.to_numeric(summary_df.get('quast_total_length'), errors='coerce')
        summary_df['cov_quast_len'] = summary_df['total_bases'] / quast_len

        summary_df['cov_quast_depth'] = pd.to_numeric(summary_df.get('quast_avg._coverage_depth'), errors='coerce')

        mash_len = pd.to_numeric(summary_df.get('mash_err_genome_size'), errors='coerce')
        summary_df['cov_mash_len'] = summary_df['total_bases'] / mash_len

        # 1. Start with Preferred: Expected Size
        summary_df['final_coverage'] = summary_df['cov_expected']
        # 2. Fallback 1: Quast Genome Size
        summary_df['final_coverage'] = summary_df['final_coverage'].fillna(summary_df['cov_quast_len'])
        # Fallback 2: Quast Avg Coverage Depth
        summary_df['final_coverage'] = summary_df['final_coverage'].fillna(summary_df['cov_quast_depth'])
        # Fallback 3: Mash Estimates
        summary_df['final_coverage'] = summary_df['final_coverage'].fillna(summary_df['cov_mash_len'])    
        summary_df['final_coverage'] = summary_df['final_coverage'].where(is_paired_end, np.nan)
        temp_cols = ['cov_expected', 'cov_quast_len', 'cov_quast_depth', 'cov_mash_len']
        summary_df = summary_df.drop(columns=[c for c in temp_cols if c in summary_df.columns])

    summary_df['coverage'] = pd.to_numeric(summary_df['final_coverage'], errors='coerce').round(2)

    return summary_df