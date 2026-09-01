import pandas as pd

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

    size_cols = [
        'checkm2_genome_size',
        'kleborate_total_size',
        'mash_err_genome_size',
        'quast_total_length',
        'expected_size'
    ]

    # Keep only the columns that actually exist in the dataframe right now
    existing_cols = [col for col in size_cols if col in summary_df.columns]

    # If we don't have at least 2 columns to compare, we can skip the check
    if len(existing_cols) >= 2:

        # Extract the sizes and safely convert to numeric (blanks/errors become NaN)
        temp_sizes = summary_df[existing_cols].apply(pd.to_numeric, errors='coerce')

        # Count how many valid, non-NaN values exist per row
        valid_counts = temp_sizes.notna().sum(axis=1)

        # Get the minimum and maximum size estimates for each row
        max_size = temp_sizes.max(axis=1)
        min_size = temp_sizes.min(axis=1)

        # Create the boolean mask: Needs at least 2 values, min_size > 0 to avoid division by zero,
        # and the difference between max and min must be greater than our threshold.
        mask_disparity = (
            (valid_counts >= 2) & 
            (min_size > 0) & 
            ((max_size - min_size) / min_size > 0.20)
        )

        # Define the warning message
        warn_msg = f"Genome size disparity (>{int(0.20 * 100)}%)"

        # Apply the warning, separated by a comma if other warnings already exist
        summary_df.loc[mask_disparity, 'warnings'] += np.where(
            summary_df.loc[mask_disparity, 'warnings'] == '',
            warn_msg,
            ',' + warn_msg
        )

    summary_df = add_quast_warnings(summary_df)

    summary_df = add_fastqc_warnings(summary_df, min_seqs=500000, max_flagged_pct=5.0, min_len=100)

    summary_df = add_checkm2_warnings(summary_df, min_completeness=90.0, max_contamination=5.0)

    summary_df = add_taxonomy_warnings(summary_df)

    summary_df = add_sylph_warnings(summary_df, min_abundance=85.0, min_ani=95.0)

    summary_df = add_st_mismatch_warning(summary_df)

    # Identify existing SeqSero note columns
    note_cols = ['seqsero_note', 'seqsero2_note', 'seqsero2s_note']
    existing_note_cols = [c for c in note_cols if c in summary_df.columns]

    if existing_note_cols:
        def has_note_text(row):
            for col in existing_note_cols:
                val = str(row[col]).strip()
                # Ignore NaNs and common empty string indicators
                if val and val.lower() not in ['nan', 'none', '-', '', '.']:
                    return True
            return False

        # Create the mask for rows containing notes
        mask_has_note = summary_df.apply(has_note_text, axis=1)

        # Define and apply the warning message
        warn_msg = "SeqSero note detected"
        
        summary_df.loc[mask_has_note, 'warnings'] = summary_df.loc[mask_has_note, 'warnings'].apply(
            lambda x: warn_msg if x == '' else f"{x},{warn_msg}"
        )


    if 'kleborate_qc' in summary_df.columns:
        mask_bad_qc = summary_df['kleborate_qc'].str.contains('fail|unreliable', case=False, na=False)
            
        warn_msg = "Kleborate QC Unreliable"
            
        summary_df.loc[mask_bad_qc, 'warnings'] = summary_df.loc[mask_bad_qc, 'warnings'].apply(
                lambda x: warn_msg if x == '' else f"{x},{warn_msg}"
            )
    return summary_df