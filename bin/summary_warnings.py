import pandas as pd
import numpy as np


def append_warning(df, mask, warn_msg, warning_column='warnings'):
    # Make sure the warning column exists
    if warning_column not in df.columns:
        df[warning_column] = ''

    # Only operate on rows where the mask is True
    mask = mask.fillna(False)

    df.loc[mask, warning_column] = np.where(
        df.loc[mask, warning_column] == '',
        warn_msg,
        df.loc[mask, warning_column] + ', ' + warn_msg
    )

    return df

def add_st_mismatch_warning(df):
    st_patterns = ['mlst_st', 'elgato_st', 'kleborate_st', 'meningotype_mlst', 'seqsero2s_st', 'shigapass_mlst']
    target_cols = [c for c in df.columns if c in st_patterns]

    if not target_cols:
        return df

    def standardize(val):
        if pd.isna(val):
            return None
        s = str(val).strip().upper()
        if s in ['', '-', 'ND', 'UNKNOWN', 'NOT FOUND', 'NONE', 'NAN', '.']:
            return None
        s = s.replace('ST', '').strip()
        if s.endswith('.0'):
            s = s[:-2]
        return s
    
    def has_mismatch(row):
        found_values = [standardize(row[col]) for col in target_cols if standardize(row[col])]
        return len(set(found_values)) > 1
    
    mask = df.apply(has_mismatch, axis=1)
    warn_msg = "ST mismatch detected"
    df.loc[mask, 'warnings'] = df.loc[mask, 'warnings'].apply(
        lambda x: warn_msg if x == '' else f"{x},{warn_msg}"
    )

    return df

def add_taxonomy_warnings(df):
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '', warn_msg, ',' + warn_msg
        )

    thresholds = {
        'kraken2_predicted_organisms': 5,
        'mash_screen_predicted_organisms': 5,
        'skani_predicted_organisms': 7,
        'sylph_predicted_organisms': 5,
        'mash_dist_predicted_organisms': 15 
    }

    for col, limit in thresholds.items():
        if col in df.columns:
            hits = pd.to_numeric(df[col], errors='coerce')
            mask_too_many = hits > limit
            tool_name = col.split('_')[0].capitalize()
            if tool_name == 'Mash':
                tool_name = "Mash " + col.split('_')[1].capitalize()
                
            warn_msg = f"Multiple organisms detected by {tool_name} (>{limit})"
            append_warning(mask_too_many, warn_msg)
    
    return df

def add seqsero2_warnings():
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


def add_kleborate_warnings():

    if 'kleborate_qc' in summary_df.columns:
        mask_bad_qc = summary_df['kleborate_qc'].str.contains('fail|unreliable', case=False, na=False)
            
        warn_msg = "Kleborate QC Unreliable"
            
        summary_df.loc[mask_bad_qc, 'warnings'] = summary_df.loc[mask_bad_qc, 'warnings'].apply(
                lambda x: warn_msg if x == '' else f"{x},{warn_msg}"
            )
    return summary_df

def add_checkm2_warnings(df, min_completeness=90.0, max_contamination=5.0):
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '', warn_msg, ',' + warn_msg
        )

    if 'checkm2_completeness' in df.columns:
        completeness = pd.to_numeric(df['checkm2_completeness'], errors='coerce')
        append_warning(completeness < min_completeness, f"Low completeness (<{min_completeness}%)")

    if 'checkm2_contamination' in df.columns:
        contamination = pd.to_numeric(df['checkm2_contamination'], errors='coerce')
        append_warning(contamination > max_contamination, f"High contamination (>{max_contamination}%)")
    
    return df

def add_genome_size_warnings(summary_df):


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


 


