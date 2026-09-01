#!/usr/bin/env python3

import numpy as np
import pandas as pd

def add_quast_warnings(df):
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '', warn_msg, ',' + warn_msg
        )

    if 'quast_n50' in df.columns:
        n50_vals = pd.to_numeric(df['quast_n50'], errors='coerce')
        append_warning(n50_vals < 30000, "Low N50 (<30000)")

    if 'quast_#_contigs' in df.columns:
        contigs = pd.to_numeric(df['quast_#_contigs'], errors='coerce')
        append_warning(contigs > 500, "High contig count (>500)")

    if 'quast_mapped_(%)' in df.columns:
        mapped = pd.to_numeric(df['quast_mapped_(%)'], errors='coerce')
        append_warning(mapped < 90.0, "Low mapping rate (<90%)")

    if "quast_#_n's_per_100_kbp" in df.columns:
        ns = pd.to_numeric(df["quast_#_n's_per_100_kbp"], errors='coerce')
        append_warning(ns > 50, "High ambiguous bases (>50 Ns/100kbp)")

    if 'quast_largest_contig' in df.columns:
        largest = pd.to_numeric(df['quast_largest_contig'], errors='coerce')
        append_warning(largest < 100000, "Largest contig too small (<100kb)")
    
    return df

def add_fastqc_warnings(df, min_seqs=500000, max_flagged_pct=5.0, min_len=100):
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '', warn_msg, ',' + warn_msg
        )

    if 'fastqc_total_sequences' in df.columns:
        total_seqs = pd.to_numeric(df['fastqc_total_sequences'], errors='coerce')
        append_warning(total_seqs < min_seqs, f"Low total sequences (<{min_seqs})")

    if 'fastqc_percent_flagged' in df.columns:
        flagged_pct = pd.to_numeric(df['fastqc_percent_flagged'], errors='coerce')
        append_warning(flagged_pct > max_flagged_pct, f"High flagged sequences (>{max_flagged_pct}%)")

    if 'fastqc_avg_length' in df.columns:
        avg_len = pd.to_numeric(df['fastqc_avg_length'], errors='coerce')
        append_warning(avg_len < min_len, f"Short avg read length (<{min_len}bp)")

    return df

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

def add_sylph_warnings(df, min_abundance=90.0, min_ani=95.0):
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '', warn_msg, ',' + warn_msg
        )

    if 'sylph_Sequence_abundance' in df.columns:
        abundance = pd.to_numeric(df['sylph_Sequence_abundance'], errors='coerce')
        append_warning(abundance < min_abundance, f"Low Sylph abundance (<{min_abundance}%)")

    if 'sylph_Adjusted_ANI' in df.columns:
        ani = pd.to_numeric(df['sylph_Adjusted_ANI'], errors='coerce')
        append_warning(ani < min_ani, f"Low Sylph ANI (<{min_ani}%)")
    
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