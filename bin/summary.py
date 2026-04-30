#!/bin/python

##########################################
# written by Erin Young                  #
# for creating summary file for grandeur #
##########################################

import pandas as pd
import json
import re
import numpy as np
import os
from os.path import exists

##########################################
# helper functions                       #
##########################################

def parse_concatenated_json(json_string):
    """Generator to parse concatenated JSON objects from a string."""
    decoder = json.JSONDecoder()
    json_string = json_string.strip()
    while json_string:
        obj, index = decoder.raw_decode(json_string)
        yield obj
        json_string = json_string[index:].lstrip()

def parse_file(summary_df, file, delim):
    print("Adding results for " + file)
    
    analysis = str(file).split("_")[0]
    new_df = pd.read_csv(file, dtype = str, index_col= False, delimiter=delim)
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = new_df.columns.str.replace('Sample', 'sample', regex=True)
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

    return summary_df


def add_quast_warnings(df):
    """
    Checks various QUAST metrics and appends warnings to the 'warnings' column:
    - Low N50 (< threshold, default 30,000)
    - High contig count (> 500)
    - Low read mapping rate (< 90%)
    - High ambiguous bases (> 50 N's per 100 kbp)
    - Largest contig too small (< 100,000 bp)
    """
    # Helper function to smoothly append warnings with a comma separator
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '',
            warn_msg,
            ',' + warn_msg
        )

    # Check: Low N50
    if 'quast_n50' in df.columns:
        n50_vals = pd.to_numeric(df['quast_n50'], errors='coerce')
        mask_low_n50 = n50_vals < 30000
        append_warning(mask_low_n50, "Low N50 (<30000)")

    # Check: High Number of Contigs (> 500)
    if 'quast_#_contigs' in df.columns:
        contigs = pd.to_numeric(df['quast_#_contigs'], errors='coerce')
        mask_high_contigs = contigs > 500
        append_warning(mask_high_contigs, "High contig count (>500)")

    # Check: Low Read Mapping Rate (< 90%)
    if 'quast_mapped_(%)' in df.columns:
        mapped = pd.to_numeric(df['quast_mapped_(%)'], errors='coerce')
        mask_low_map = mapped < 90.0
        append_warning(mask_low_map, "Low mapping rate (<90%)")

    # Check: High Ambiguous Bases (> 50 N's per 100 kbp)
    if "quast_#_n's_per_100_kbp" in df.columns:
        ns = pd.to_numeric(df["quast_#_n's_per_100_kbp"], errors='coerce')
        mask_high_ns = ns > 50
        append_warning(mask_high_ns, "High ambiguous bases (>50 Ns/100kbp)")

    # Check: Small Largest Contig (< 100,000 bp)
    if 'quast_largest_contig' in df.columns:
        largest = pd.to_numeric(df['quast_largest_contig'], errors='coerce')
        mask_small_largest = largest < 100000
        append_warning(mask_small_largest, "Largest contig too small (<100kb)")
    
    return df

def add_fastqc_warnings(df, min_seqs=500000, max_flagged_pct=5.0, min_len=100):
    """
    Checks FastQC metrics and appends warnings to the 'warnings' column:
    - Low total sequences (< min_seqs, default 500,000)
    - High percent flagged sequences (> max_flagged_pct, default 5.0%)
    - Short average read length (< min_len, default 100 bp)
    """
    # Helper function to smoothly append warnings with a comma separator
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '',
            warn_msg,
            ',' + warn_msg
        )

    # Check: Low Total Sequences
    if 'fastqc_total_sequences' in df.columns:
        total_seqs = pd.to_numeric(df['fastqc_total_sequences'], errors='coerce')
        mask_low_seqs = total_seqs < min_seqs
        append_warning(mask_low_seqs, f"Low total sequences (<{min_seqs})")

    # Check: High Flagged Sequences
    if 'fastqc_percent_flagged' in df.columns:
        flagged_pct = pd.to_numeric(df['fastqc_percent_flagged'], errors='coerce')
        mask_high_flagged = flagged_pct > max_flagged_pct
        append_warning(mask_high_flagged, f"High flagged sequences (>{max_flagged_pct}%)")

    # Check: Short Average Length
    if 'fastqc_avg_length' in df.columns:
        avg_len = pd.to_numeric(df['fastqc_avg_length'], errors='coerce')
        mask_short_len = avg_len < min_len
        append_warning(mask_short_len, f"Short avg read length (<{min_len}bp)")

    return df

def add_checkm2_warnings(df, min_completeness=90.0, max_contamination=5.0):
    """
    Checks CheckM2 metrics and appends warnings to the 'warnings' column:
    - Low completeness (< min_completeness, default 90.0%)
    - High contamination (> max_contamination, default 5.0%)
    """
    # Helper function to smoothly append warnings with a comma separator
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '',
            warn_msg,
            ',' + warn_msg
        )

    # Check: Low Completeness
    if 'checkm2_completeness' in df.columns:
        completeness = pd.to_numeric(df['checkm2_completeness'], errors='coerce')
        mask_low_comp = completeness < min_completeness
        append_warning(mask_low_comp, f"Low completeness (<{min_completeness}%)")

    # Check: High Contamination
    if 'checkm2_contamination' in df.columns:
        contamination = pd.to_numeric(df['checkm2_contamination'], errors='coerce')
        mask_high_contam = contamination > max_contamination
        append_warning(mask_high_contam, f"High contamination (>{max_contamination}%)")
    
    return df

def add_taxonomy_warnings(df):
    """
    Checks taxonomic classification columns and appends a warning if the 
    number of predicted organisms exceeds the tool's specific threshold.
    """

    # Helper function to smoothly append warnings with a comma separator
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '',
            warn_msg,
            ',' + warn_msg
        )

    # A threshold of 1 means anything 2 or higher gets flagged.
    thresholds = {
        'kraken2_predicted_organisms': 5,
        'mash_screen_predicted_organisms': 5,
        'skani_predicted_organisms': 7,
        'sylph_predicted_organisms': 5,
        # mash_dist tends to have higher hit counts natively
        'mash_dist_predicted_organisms': 15 
    }

    # Loop through the columns and apply the warnings
    for col, limit in thresholds.items():
        if col in df.columns:
            # Convert to numeric safely
            hits = pd.to_numeric(df[col], errors='coerce')
            
            # Mask where hits are greater than the allowed limit
            mask_too_many = hits > limit
            
            # Create a clean, tool-specific warning message
            tool_name = col.split('_')[0].capitalize() # e.g., 'Kraken2', 'Skani'
            if tool_name == 'Mash':
                tool_name = "Mash " + col.split('_')[1].capitalize() # 'Mash Dist' or 'Mash Screen'
                
            warn_msg = f"Multiple organisms detected by {tool_name} (>{limit})"
            
            append_warning(mask_too_many, warn_msg)
    
    return df

def add_sylph_warnings(df, min_abundance=90.0, min_ani=95.0):
    """
    Checks Sylph metrics and appends warnings to the 'warnings' column:
    - Low Sequence Abundance (< 90%)
    - Low Adjusted ANI (< 95%)
    """

    # Helper function to smoothly append warnings with a comma separator
    def append_warning(mask, warn_msg):
        df.loc[mask, 'warnings'] += np.where(
            df.loc[mask, 'warnings'] == '',
            warn_msg,
            ',' + warn_msg
        )

    # Check: Low Abundance (Using Sequence_abundance as it's the raw data metric)
    if 'sylph_Sequence_abundance' in df.columns:
        abundance = pd.to_numeric(df['sylph_Sequence_abundance'], errors='coerce')
        mask_low_abund = abundance < min_abundance
        append_warning(mask_low_abund, f"Low Sylph abundance (<{min_abundance}%)")

    # Check: Low Adjusted ANI
    if 'sylph_Adjusted_ANI' in df.columns:
        ani = pd.to_numeric(df['sylph_Adjusted_ANI'], errors='coerce')
        mask_low_ani = ani < min_ani
        append_warning(mask_low_ani, f"Low Sylph ANI (<{min_ani}%)")
    
    return df


def add_st_mismatch_warning(df):
    """
    Compares all available ST columns and adds a warning if 
    conflicting STs are found within a single row.
    """

    st_patterns = ['mlst_st', 'elgato_st', 'kleborate_st', 'meningotype_mlst', 'seqsero2s_st', 'shigapass_mlst']
    
    target_cols = [c for c in df.columns if c in st_patterns]

    if not target_cols:
        return df

    def standardize(val):
        """Cleans ST values for direct comparison."""
        if pd.isna(val):
            return None
        s = str(val).strip().upper()
        # Filter out common empty/null placeholders
        if s in ['', '-', 'ND', 'UNKNOWN', 'NOT FOUND', 'NONE', 'NAN', '.']:
            return None
        
        # Remove 'ST' prefix and any trailing '.0' from floats
        s = s.replace('ST', '').strip()
        if s.endswith('.0'):
            s = s[:-2]
        return s
    
    def has_mismatch(row):
        # Collect all unique, cleaned ST values present in this row
        found_values = []
        for col in target_cols:
            clean = standardize(row[col])
            if clean:
                found_values.append(clean)
        
        # If we have multiple unique values (e.g., {'11', '14464'}), it's a mismatch
        return len(set(found_values)) > 1
    
    # 3. Identify rows with conflicts
    mask = df.apply(has_mismatch, axis=1)

    # 4. Append warning with comma handling
    warn_msg = "ST mismatch detected"
    df.loc[mask, 'warnings'] = df.loc[mask, 'warnings'].apply(
        lambda x: warn_msg if x == '' else f"{x},{warn_msg}"
    )

    return df


def create_node(parent=None):
    """Creates a basic dictionary representing a node."""
    return {
        'parent': parent,
        'children': [],
        'name': "",
        'length': 0.0
    }

def parse_newick(newick_str):
    """Parses a newick string into linked dictionaries."""
    newick_str = newick_str.strip().replace('\n', '').replace('\r', '')
    root = create_node()
    current_node = root
    state = 'name'
    buffer = ""

    for char in newick_str:
        if char == '(':
            new_node = create_node(parent=current_node)
            current_node['children'].append(new_node)
            current_node = new_node
            state = 'name'
            buffer = ""
        elif char == ',':
            if state == 'name':
                current_node['name'] = buffer.strip()
            else:
                if buffer.strip():
                    current_node['length'] = float(buffer)
            
            new_node = create_node(parent=current_node['parent'])
            current_node['parent']['children'].append(new_node)
            current_node = new_node
            state = 'name'
            buffer = ""
        elif char == ')':
            if state == 'name':
                current_node['name'] = buffer.strip()
            else:
                if buffer.strip():
                    current_node['length'] = float(buffer)
            
            current_node = current_node['parent']
            state = 'name'
            buffer = ""
        elif char == ':':
            if state == 'name':
                current_node['name'] = buffer.strip()
            state = 'length'
            buffer = ""
        elif char == ';':
            if state == 'name':
                current_node['name'] = buffer.strip()
            elif state == 'length':
                if buffer.strip():
                    current_node['length'] = float(buffer)
            break
        else:
            buffer += char

    return root

def get_tip_distance_stats(newick):
    """Calculates the average, min, and max distances using dictionaries."""

    with open(newick, 'r') as f:
        newick_str = f.read()

    file_basename = os.path.basename(newick)

    root = parse_newick(newick_str)
    
    tips = []
    def find_tips(node):
        if len(node['children']) == 0 and node['name']: 
            tips.append(node)
        for child in node['children']:
            find_tips(child)
            
    find_tips(root)
    
    if len(tips) < 2:
        return pd.DataFrame()
        
    results = {}
    
    for start_tip in tips:
        distances = {}
        visited = set()
        queue = [(start_tip, 0.0)]
        visited.add(id(start_tip))
        
        while queue:
            current, dist = queue.pop(0)
            is_tip = len(current['children']) == 0
            
            if is_tip and id(current) != id(start_tip):
                distances[current['name']] = dist
                
            parent = current['parent']
            if parent and id(parent) not in visited:
                visited.add(id(parent))
                queue.append((parent, dist + current['length']))
                
            for child in current['children']:
                if id(child) not in visited:
                    visited.add(id(child))
                    queue.append((child, dist + child['length']))
                    
        dist_values = list(distances.values())
        results[start_tip['name']] = {
            'average': sum(dist_values) / len(dist_values),
            'min': min(dist_values),
            'max': max(dist_values)
        }
    
    # 3. Add the basename to the dictionary comprehensions
    df_data = [{'file': file_basename, 'sample': tip, **metrics} for tip, metrics in results.items()]
    
    # 4. Convert to a Pandas DataFrame
    df = pd.DataFrame(df_data)
    df = df.add_prefix(file_basename + "_")
    df['sample'] = df[file_basename + "_sample"]

    return df


def get_snp_distance_stats(filepath):
    """Calculates the average, min, and max distances from a snp-dists matrix."""
    
    # 1. Read the matrix. 
    # index_col=0 tells pandas to use the first column (sample names) as the row labels.
    df_matrix = pd.read_csv(filepath, index_col=0)
    
    # 2. Extract the basename of the file
    file_basename = os.path.basename(filepath)
    
    results = []
    
    # 3. Loop through each sample in the index
    for sample in df_matrix.index:
        # Get the distances for this sample, and drop the distance to itself
        # We ensure 'sample' is treated as a string to match pandas column names
        distances = df_matrix.loc[sample].drop(str(sample))
        
        # Calculate the stats and append to our results list
        results.append({
            'file': file_basename,
            'sample': sample,
            'average': distances.mean(),
            'min': distances.min(),
            'max': distances.max()
        })
        
    # 4. Convert to a Pandas DataFrame
    df = pd.DataFrame(results)
    df = df.add_prefix(file_basename + "_")
    df['sample'] = df[file_basename + "_sample"]
    
    return df

##########################################
# defining files                         #
##########################################

# input files
names          = 'input_files.txt'
amrfinderplus  = 'amrfinderplus_summary.txt'
checkm2        = 'checkm2_summary.tsv'
core           = 'multiqc_core_genome_evaluation.txt'
datasets       = 'datasets_summary.csv'
drprg          = 'drprg_summary.tsv'
elgato         = 'elgato_summary.tsv'
emmtyper       = 'emmtyper_summary.tsv'
skani          = 'skani_summary.tsv'
fastqc         = 'fastqc_summary.csv'
genome_sizes   = "genome_sizes.json"
kaptive        = "kaptive_summary.tsv"
kleborate      = 'kleborate_results.tsv'
kraken2        = 'kraken2_summary.csv'
legsta         = 'legsta_summary.csv'
mash_dist      = 'mashdist_summary.csv'
mash_screen    = 'mashscreen_summary.txt'
mash_err       = 'mash_err_summary.csv'
meningotype    = 'meningotype_summary.tsv'
mlst           = 'mlst_summary.tsv'
mykrobe        = 'mykrobe_summary.csv'
ngmaster       = 'ngmaster_summary.csv'
pbptyper       = 'pbptyper_summary.tsv'
plasmidfinder  = 'plasmidfinder_result.json'
quast          = 'quast_report.tsv'
quast_contig   = 'quast_contig_report.tsv'
seqsero2       = 'seqsero2_results.txt'
seqsero2s      = 'seqsero2s_results.txt'
serotypefinder = 'serotypefinder_results.txt'
shigapass      = 'shigapass_summary.tsv'
spestimator    = 'spestimator_summary.csv'
sylph          = 'sylph_summary.tsv'
multiqc_json   = 'multiqc_data.json'
multiqc_stats  = 'multiqc_general_stats.txt'

gotree = 'gotree_summary.tsv'
snpdist_matrices = ['snpdists_core_gene_alignment.txt', 'snpdists_ska_alignment.txt']
newick_files = ['iqtree_core_gene_alignment.treefile.nwk', 'iqtree_ska_alignment.treefile.nwk', 'mashtree.nwk']

# final files
final          = 'grandeur_summary'
extended       = 'summary/grandeur_extended_summary'

##########################################
# grouping similar files                 #
##########################################

csv_files = [ legsta, mykrobe, ngmaster ]
tsv_files = [ drprg, checkm2, elgato, meningotype, seqsero2, seqsero2s, shigapass, kaptive, kleborate, mlst, emmtyper, pbptyper ]

##########################################
# exiting if no input files              #
##########################################

if not exists(names) :
    print("No analyses to report on for this run!")
    with open(extended + '.tsv', 'w') as fp:
        pass
    with open(extended + '.txt', 'w') as fp:
        pass
    with open(final + '.tsv', 'w') as fp:
        pass
    with open(final + '.txt', 'w') as fp:
        pass
    quit()

##########################################
# creating the summary dataframe         #
##########################################

input_cols = ['sample', 'file', 'file_2', 'version']

summary_df = pd.read_csv(names, dtype = str, names=input_cols, delimiter=",")
summary_df['sample'] = summary_df['sample'].astype(str)
summary_df['file'] = summary_df['file'].astype(str)
summary_df['warnings'] = ''
columns = list(summary_df.columns)

# csv files
for file in csv_files :
    if exists(file) :
        summary_df = parse_file(summary_df, file, ",")

# tsv files
for file in tsv_files :
    if exists(file) :
        summary_df = parse_file(summary_df, file, "\t")        

# for specific tools

# amrfinderplus : merging many rows into one with relevant information
if exists(amrfinderplus):
    print("Adding results for " + amrfinderplus)
    analysis = "amrfinder"
    amr_df = pd.read_table(amrfinderplus, dtype=str, index_col=False)

    # --- 1. PREP FOR SUMMARY_DF (List of genes) ---
    summary_prep = amr_df.copy()
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
    summary_df = pd.merge(summary_df, summary_genes, on="sample", how='left')

    # --- 3. STANDALONE MATRIX FILE ---
    # Create the "Cov,Ident" value string for the matrix
    amr_df['matrix_val'] = amr_df['% Coverage of reference'] + ',' + amr_df['% Identity to reference']

    # Create the detailed header "Type_Family_Gene"
    def format_gene_header(row):
        if 'family' in str(row['Element name']).lower():
            family = str(row['Element name']).split(' family')[0] + " family"
            return f"{row['Type']}_{family}_{row['Element symbol']}"
        return f"{row['Type']}_{row['Element symbol']}"

    amr_df['gene_header'] = amr_df.apply(format_gene_header, axis=1)
    # TODO : finish creating AMR family csv file
    # Pivot to create the wide matrix
    # amr_matrix = amr_df.pivot(index='Name', columns='gene_header', values='matrix_val').fillna('-')
    
    # Save the matrix to a separate file
    # amr_matrix.to_csv('amr_virulence_matrix.tsv', sep='\t')
    # print("Standalone AMR matrix saved to amr_virulence_matrix.tsv")

# datasets : adding a count of reference genomes available
if exists(datasets):
    print("Adding reference genome counts from " + datasets)
    # Read the datasets file
    datasets_df = pd.read_csv(datasets, dtype=str)
    
    # Calculate the total number of genomes in the database
    num_genomes = len(datasets_df)
    
    # Assign this constant value to a new column in the main summary
    summary_df['datasets_num_genomes'] = num_genomes


# (Assuming summary_df and fastqc variables are defined)
if exists(fastqc):
    file = fastqc
    print("Adding results for " + file)
    analysis = "fastqc"
    new_df = pd.read_csv(file, dtype=str, index_col=False)
    
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

# kraken2 : merging relevant rows into one
if exists(kraken2):
    print("Adding results for " + kraken2)
    analysis = "kraken2"
    new_df = pd.read_csv(kraken2, dtype=str, index_col=False)
    
    # Sort by abundance to find the top hit
    new_df = new_df.sort_values(['Sample', 'Percentage of fragments'], ascending=False)
    
    # Save the top organism name for the "Predicted Organism" logic later
    tmp_df = new_df.drop_duplicates(subset=['Sample'], keep="first")[['Sample', 'Scientific name']]
    tmp_df.columns = ['Sample', 'kraken2_top_organism']

    # Create the string format for the report (e.g., "E.coli (98%)")
    new_df['kraken2_formatted'] = new_df['Scientific name'] + " (" + new_df['Percentage of fragments'] + "%)"
    
    # Group by Sample and calculate both the report list and the unique organism count
    grouped = new_df.groupby('Sample').agg(
        kraken2_report=('kraken2_formatted', lambda x: ', '.join(x.dropna().astype(str))),
        kraken2_predicted_organisms=('Scientific name', 'nunique')
    ).reset_index()
    
    # Merge both dataframes into the main summary_df
    summary_df = pd.merge(summary_df, tmp_df, left_on="sample", right_on="Sample", how='left').drop('Sample', axis=1)
    summary_df = pd.merge(summary_df, grouped, left_on="sample", right_on="Sample", how='left').drop('Sample', axis=1)

# mash dist
if exists(mash_dist) :
    file = mash_dist
    print("Adding results for " + file)
    analysis = "mash_dist"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    # header : sample,reference,query,mash-distance,P-value,matching-hashes,organism
    new_df = new_df.sort_values(by = ['sample', 'P-value', 'mash-distance'], ascending = [True, True, True])

    counts = new_df.groupby('sample')['organism'].nunique().reset_index()
    counts.columns = ['sample', 'predicted_organisms']

    new_df = pd.merge(new_df, counts, on='sample', how='left')
    new_df = new_df.drop_duplicates(subset=['sample'], keep='first').reset_index(drop=True)
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# mash screen
if exists(mash_screen) :
    file = mash_screen
    print("Adding results for " + file)
    analysis = "mash_screen"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    # header : sample,identity,shared-hashes,median-multiplicity,p-value,query-ID,organism
    new_df = new_df.sort_values(by = ['sample', 'p-value', 'identity'], ascending = [True, True, False])

    counts = new_df.groupby('sample')['organism'].nunique().reset_index()
    counts.columns = ['sample', 'predicted_organisms']

    new_df = new_df.drop_duplicates(subset=['sample'], keep='first')
    new_df = pd.merge(new_df, counts, on='sample', how='left')
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

if exists(mash_err):
    file = mash_err
    print("Adding results for " + file)
    analysis = "mash_err"
    
    with open(file, 'r') as f:
        content = f.read()
    
    pattern = re.compile(
        r"Estimated genome size:\s+([\d.e+]+).*?"
        r"Estimated coverage:\s+([\d.]+).*?"
        r"Writing to\s+(.*?)\.msh", 
        re.DOTALL
    )
    
    matches = pattern.findall(content)
    
    data = []
    for size, coverage, sample_id in matches:
        data.append({
            "sample": sample_id.strip(),
            "genome_size": float(size),
            "coverage": float(coverage)
        })
    
    # Create DataFrame
    new_df = pd.DataFrame(data)
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# plasmidfinder : merging relevant rows into one
if exists(plasmidfinder) :
    file = plasmidfinder
    print("Adding results for " + file)
    analysis = "plasmidfinder"
    with open(file, 'r') as f:
        content = f.read()

    all_data = []
    
    for run in parse_concatenated_json(content):
        # Extract sample name from the JSON value
        execs = run.get("software_executions", {})
        first_exec = next(iter(execs.values())) if execs else {}
        out_json_val = first_exec.get("parameters", {}).get("out_json", "")
        sample_name = out_json_val.split("/")[-2]

        seq_regions = run.get("seq_regions", {})
        
        plasmids = ", ".join([f"{d.get('name', '')} ({d.get('identity', '')})" for d in seq_regions.values()])
        all_data.append({
            "sample": sample_name,
            "plasmid_(identity)": plasmids if plasmids else ""
        })

    new_df = pd.DataFrame(all_data)
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = new_df.columns.str.replace('Sample', 'sample', regex=True)
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# quast : combining both files
q_df  = pd.DataFrame()
qc_df = pd.DataFrame()
if exists(quast):
    print("Adding results for " + quast)
    file = quast
    analysis = str(file).split("_")[0]
    q_df = pd.read_table(file, dtype = str, index_col= False)
    q_df = q_df.add_prefix(analysis + "_")
    q_df.columns = [x.lower() for x in q_df.columns]

if exists(quast_contig):
    print("Adding results for " + quast_contig)
    file = quast_contig
    analysis = str(file).split("_")[0]
    qc_df = pd.read_table(file, dtype = str, index_col= False)
    qc_df = qc_df.add_prefix(analysis + "_")
    qc_df.columns = [x.lower() for x in qc_df.columns]

if exists(quast) or exists(quast_contig):
    new_df = pd.concat([q_df, qc_df])

    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how='left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# serotypefinder : splitting O and H groups, getting the top hit for O and H group, combining rows
if exists(serotypefinder):
    file = serotypefinder
    print("Adding results for " + file)
    analysis = "serotypefinder"
    new_df = pd.read_table(file, dtype=str, index_col=False)
    
    counts_df = new_df.groupby(['sample', 'Database']).size().unstack(fill_value=0).reset_index()
    counts_df = counts_df.rename_axis(None, axis=1) # Clean up the column grouping name
    
    if 'O_type' not in counts_df.columns:
        counts_df['O_type'] = 0
    if 'H_type' not in counts_df.columns:
        counts_df['H_type'] = 0
        
    counts_df = counts_df.rename(columns={
        'O_type': analysis + '_O_count', 
        'H_type': analysis + '_H_count'
    })

    new_df = new_df.sort_values(by='Identity', ascending=False)
    new_df = new_df.drop_duplicates(subset=['sample', 'Database'], keep="first")
    new_df = new_df.add_prefix(analysis + '_')
    
    H_df = new_df[new_df[analysis + '_Database' ] == 'H_type'].copy()
    H_df = H_df.add_suffix('_H')
    
    O_df = new_df[new_df[analysis + '_Database' ] == 'O_type'].copy()
    O_df = O_df.add_suffix('_O')
    
    summary_df = pd.merge(summary_df, O_df, left_on="sample", right_on=analysis + "_sample_O", how='left')
    summary_df.drop(analysis + "_sample_O", axis=1, inplace=True)
    
    summary_df = pd.merge(summary_df, H_df, left_on="sample", right_on=analysis + "_sample_H", how='left')
    summary_df.drop(analysis + "_sample_H", axis=1, inplace=True)
    
    summary_df = pd.merge(summary_df, counts_df[['sample', analysis + '_O_count', analysis + '_H_count']], on="sample", how="left")


# skani
if exists(skani):
    print("Adding results for " + skani)
    analysis = "skani"
    new_df = pd.read_table(skani, dtype=str, index_col=False)
    
    # Clean up names and split Genus_species
    new_df['organism'] = new_df['Ref_file'].str.split('_').str[0:2].str.join('_')
    new_df = new_df.sort_values(['sample', 'ANI'], ascending=[True, False])
    new_df = new_df.add_prefix('skani_')

    top_df = new_df.drop_duplicates(subset=['skani_sample'], keep='first').copy()
    
    # Count how many unique organisms Skani thinks are in this one "isolate"
    counts_df = new_df.groupby('skani_sample').agg(
        skani_predicted_organisms=('skani_organism', 'nunique')
    ).reset_index()
    summary_df = pd.merge(summary_df, top_df, left_on="sample", right_on="skani_sample", how='left').drop('skani_sample', axis=1)
    summary_df = pd.merge(summary_df, counts_df, left_on="sample", right_on="skani_sample", how='left')

# spestimator : counting unique reference hits per sample
if exists(spestimator):
    print("Adding results for " + spestimator)
    analysis = "spestimator"
    
    # Read the TSV file
    new_df = pd.read_table(spestimator, sep=",", dtype=str)
    
    # 1. Clean the 'input file' column to match your 'sample' IDs
    # This removes '_contigs.fa' and other extensions
    new_df['sample'] = (
        new_df['input file']
        .str.replace('_contigs.fa', '', regex=False)
        .str.replace(r'\.(fasta|fna|fa)$', '', regex=True)
    )
    
    # 2. Count unique organisms identified for each sample
    # Using 'organism' as the unique identifier here
    sp_counts = new_df.groupby('sample')['organism'].nunique().reset_index()
    sp_counts.columns = ['sample', 'spestimator_num_refs']
    
    # 3. Merge into the main summary_df
    summary_df = pd.merge(summary_df, sp_counts, on='sample', how='left')

# sylph
if exists(sylph):
    print("Adding results for " + sylph)
    new_df = pd.read_table(sylph, dtype=str, index_col=False)
    
    # Convert abundance to numeric so we can safely sort by it
    new_df['Taxonomic_abundance'] = pd.to_numeric(new_df['Taxonomic_abundance'], errors='coerce')
    
    # Sort to ensure the top hit (highest abundance) is first for every sample
    new_df = new_df.sort_values(['sample', 'Taxonomic_abundance'], ascending=[True, False])
    
    # --- 1. Get the Top Hit ---
    # Keep only the first row per sample to maintain the main dataframe columns
    top_df = new_df.drop_duplicates(subset=['sample'], keep='first').copy()
    top_df = top_df.add_prefix('sylph_')
    
    # --- 2. Count Unique Organisms ---
    # Group by sample and count the unique Genome_files
    counts_df = new_df.groupby('sample').agg(
        sylph_predicted_organisms=('Genome_file', 'nunique')
    ).reset_index()
    
    # --- 3. Merge ---
    # Merge the top hit details
    summary_df = pd.merge(summary_df, top_df, left_on="sample", right_on="sylph_sample", how='left').drop('sylph_sample', axis=1)
    
    # Merge the organism count
    summary_df = pd.merge(summary_df, counts_df, on="sample", how='left')

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

# core genome analysis file is also from multiqc
if exists(core):
    file = core
    analysis = "core_genome_genes"
    print("Adding core genome percentage from " + file)
    new_df = pd.read_table(file, dtype = str, index_col= False)
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_Sample", how = 'left')
    summary_df.drop(analysis + "_Sample", axis=1, inplace=True)
    summary_df['per_core_genome_genes'] = summary_df[analysis + '_core'].astype(float) / (summary_df[analysis + '_soft'].astype(float) + summary_df[analysis + '_core'].astype(float) + summary_df[analysis + '_shell'].astype(float) + summary_df[analysis + '_cloud'].astype(float))
    summary_df['per_core_genome_genes'] = summary_df['per_core_genome_genes'].astype(float) * 100
    summary_df['per_core_genome_genes'] = summary_df['per_core_genome_genes'].round(2)

##########################################
# predicting organism                    #
##########################################

print("Predicting organism")

summary_df['predicted_organism'] = pd.NA

if 'skani_organism' in summary_df.columns:
    summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['skani_organism'])

if 'kraken2_top_organism' in summary_df.columns:
    summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['kraken2_top_organism'])

if 'mash_screen_organism' in summary_df.columns:
    summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['mash_screen_organism'])

if 'mash_dist_organism' in summary_df.columns:
    summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['mash_dist_organism'])

summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna("Unknown")

summary_df['predicted_organism'] = summary_df['predicted_organism'].str.strip()

# Adjusting shigellas
if 'shigapass_predicted_serotype' in summary_df.columns:
    print("Refining E. coli and Shigella predictions using ShigaPass")
    
    def refine_shigella_ecoli(row):
        org = str(row.get('predicted_organism', ''))
        shiga_pred = str(row.get('shigapass_predicted_serotype', ''))
        
        if 'escherichia' in org.lower() or 'shigella' in org.lower():
            
            if pd.notna(shiga_pred) and shiga_pred.strip() != '' and shiga_pred.lower() != 'nan':
                shiga_pred_clean = shiga_pred.strip().upper()
                new_org = None
                
                if shiga_pred_clean.startswith('SS'):
                    new_org = 'Shigella sonnei'
                elif shiga_pred_clean.startswith('SF'):
                    new_org = 'Shigella flexneri'
                elif shiga_pred_clean.startswith('SB'):
                    new_org = 'Shigella boydii'
                elif shiga_pred_clean.startswith('SD'):
                    new_org = 'Shigella dysenteriae'
                elif 'NOT SHIGELLA' in shiga_pred_clean or 'EIEC' in shiga_pred_clean:
                    if 'shigella' in org.lower():
                        new_org = 'Escherichia coli'
                    else:
                        return org
                
                if new_org:
                    if '_' in org:
                        return new_org.replace(' ', '_')
                    else:
                        return new_org
                        
        return org

    summary_df['predicted_organism'] = summary_df.apply(refine_shigella_ecoli, axis=1)

##########################################
# size and coverage estimates            #
##########################################

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

##########################################
# summarizing phylogenetics              #
##########################################

for nwk in newick_files:
    if exists(nwk) :
        print("Adding results for " + nwk)
        new_df = get_tip_distance_stats(nwk)
        summary_df = pd.merge(summary_df, new_df, on="sample", how = 'left')

for snp_matrix in snpdist_matrices:
    if exists(snp_matrix) :
        print("Adding results for " + snp_matrix)
        new_df = get_snp_distance_stats(snp_matrix)
        summary_df = pd.merge(summary_df, new_df, on="sample", how = 'left')

if exists(gotree):
    print("Adding average values from gotree")
    new_df = pd.read_csv(gotree, sep='\t')
    for index, row in new_df.iterrows():
        analysis_name = row['sample']
        
        prefix = analysis_name.replace('gotree_', '').replace('.treefile', '')
        
        # 3. Create the new column names
        mean_col = f"{prefix}_meanbrlen"
        sum_col = f"{prefix}_sumbrlen"
        
        # 4. Assign the values to the main DataFrame
        # Pandas will automatically broadcast this single value to every row
        summary_df[mean_col] = row['meanbrlen']
        summary_df[sum_col]  = row['sumbrlen']

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


##########################################
# creating files                         #
##########################################

print("Creating final files")

summary_df = summary_df.sort_values(by='sample')
summary_df = summary_df.fillna("")

summary_df.columns = summary_df.columns.str.replace(' ', '_')

summary_df.to_csv(extended + '.tsv', index=False, sep="\t")
summary_df.to_csv(extended + '.txt', index=False, sep=";")

transposed_df = summary_df.set_index('sample').T
transposed_df = transposed_df.reset_index().rename(columns={'index': 'metric'})
transposed_df.to_csv(extended + '_transposed.tsv', index=False, sep="\t")
transposed_df.to_csv(extended + '_transposed.txt', index=False, sep=";")


# reducing to the top 1 or 2 results for each analysis
final_columns = [
    # general information
    'coverage',
    'per_core_genome_genes',
    'fastqc_total_sequences',
    'fastqc_flagged_sequences',
    'fastqc_avg_length',
    'fastp_pct_passed_reads',
    'quast_#_contigs',
    'quast_mapped_(%)',
    'quast_gc_(%)',
    'checkm2_completeness',
    'checkm2_contamination',
    'warnings',
    'amrfinder_genes_(per_cov/per_ident)',

    # species
    'predicted_organism',
    'mlst_matching_pubmlst_scheme',
    'mlst_st',
    'skani_ANI',
    'skani_organism',
    'spestimator_num_refs',
    'datasets_num_genomes',
    'sylph_Adjusted_ANI',
    'sylph_Taxonomic_abundance',
    'mash_dist_organism',
    'mash_dist_mash-distance',
    'mash_dist_p-value',
    'mash_screen_organism',
    'mash_screen_identity',
    'mash_screen_shared-hashes',
    'mash_screen_p-value',
    'plasmidfinder_plasmid_(identity)',

    # contamination
    'kraken2_report',

    # species specific information
    'seqsero2_predicted_antigenic_profile',
    'seqsero2_predicted_serotype',
    'emmtyper_predicted_emm-type',
    'kleborate_virulence_score',
    'kleborate_resistance_score',
    'kaptive_acinetobacter_baumannii_oc_locus_primary_reference_best_match_locus',
    'kaptive_acinetobacter_baumannii_k_locus_primary_reference_best_match_locus',
    'kaptive_klebsiella_k_locus_primary_reference_best_match_locus',
    'kaptive_klebsiella_o_locus_primary_reference_best_match_locus',
    'kaptive_vibriopara_kaptivedb_k_best_match_locus',
    'kaptive_vibriopara_kaptivedb_o_best_match_locus',
    'elgato_st',
    'meningotype_serogroup',
    'ngmaster_ng-mast/ng-star',
    'mykrobe_phylo_group',
    'mykrobe_species',
    'mykrobe_lineage',
    'drprg_susceptibility',
    'pbptyper_pbptype',
    'serotypefinder_Serotype_O',
    'serotypefinder_Serotype_H',
    'shigapass_ipah',
    'shigapass_predicted_serotype',

    # phylogenetic analysis results
    'snpdists_core_gene_alignment.txt_average',
    'snpdists_core_gene_alignment.txt_min',
    'snpdists_core_gene_alignment.txt_max',
    'iqtree_core_gene_alignment.treefile.nwk_average',
    'iqtree_core_gene_alignment.treefile.nwk_min',
    'iqtree_core_gene_alignment.treefile.nwk_max',
    'iqtree_core_gene_alignment_meanbrlen',
    'iqtree_core_gene_alignment_sumbrlen',
    'snpdists_ska_alignment.txt_average',
    'snpdists_ska_alignment.txt_min',
    'snpdists_ska_alignment.txt_max',
    'iqtree_ska_alignment.treefile.nwk_average',
    'iqtree_ska_alignment.treefile.nwk_min',
    'iqtree_ska_alignment.treefile.nwk_max',
    'iqtree_ska_alignment_meanbrlen',
    'iqtree_ska_alignment_sumbrlen',
    'mashtree.nwk_average',
    'mashtree.nwk_min',
    'mashtree.nwk_max',
    'mashtree_meanbrlen',
    'mashtree_sumbrlen'

    ]

set_columns = []
for new_column in final_columns :
    if new_column in summary_df.columns :
        set_columns.append(new_column)

summary_df.to_csv(final + '.tsv', columns = ['sample','file','version'] + set_columns, index=False, sep="\t")
summary_df.to_csv(final + '.txt', columns = ['sample','file','version'] + set_columns, index=False, sep=";")
