#!/bin/python

##########################################
# written by Erin Young                  #
# for creating summary file for grandeur #
##########################################

import pandas as pd
import json
import re
from os.path import exists

##########################################
# helper functions                       #
##########################################

def add_warning(df, new_warn_col):
    """Safely appends new warnings to the main 'warnings' column."""
    if new_warn_col in df.columns:
        # Fill NaNs so string concatenation doesn't result in NaN
        df['warnings'] = df['warnings'].fillna('')
        df[new_warn_col] = df[new_warn_col].fillna('')
        
        # Join with a comma only if both strings have content
        df['warnings'] = df.apply(
            lambda row: ", ".join(filter(None, [str(row['warnings']), str(row[new_warn_col])])), 
            axis=1
        )
        # Drop the tool-specific warning column to keep the dataframe clean
        df.drop(new_warn_col, axis=1, inplace=True)
    return df

def check_isolate_purity(group):
    """Checks a group of taxonomic hits for contamination signs."""
    notes = []
    # Check secondary hits for contamination (>1% abundance)
    secondary_hits = group.iloc[1:] 
    significant_contam = secondary_hits[secondary_hits['Taxonomic_abundance'] > 1.0]
    
    if not significant_contam.empty:
        notes.append(f"sylph detected {len(significant_contam)} secondary species >1%")
        
    # Check primary hit for general purity
    if not group.empty and group.iloc[0]['Taxonomic_abundance'] < 98.0:
        notes.append("sylph predicts low purity isolate (<98%)")
        
    return ", ".join(notes)



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
    new_df = new_df.replace('Sample', 'sample', regex=True)
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

    return summary_df


##########################################
# defining files                         #
##########################################

# input files
names          = 'input_files.txt'
amrfinderplus  = 'amrfinderplus.txt'
checkm2        = 'checkm2_summary.tsv'
core           = 'multiqc_core_genome_evaluation-plot.txt'
datasets       = 'datasets_summary.csv'
drprg          = 'drprg_summary.tsv'
elgato         = 'elgato_summary.tsv'
emmtyper       = 'emmtyper_summary.tsv'
skani          = 'skani_summary.tsv'
fastqc         = 'fastqc_summary.csv'
genome_sizes   = "genome_sizes.json"
kaptive        = "kaptive_summary.txt"
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
spestimator    = 'spestimator_summary.tsv'
sylph          = 'sylph_summary.tsv'
multiqc_json   = 'multiqc_data.json'
multiqc_stats  = 'multiqc_general_stats.txt'

# final files
final          = 'grandeur_summary'
extended       = 'summary/grandeur_extended_summary'

##########################################
# grouping similar files                 #
##########################################

csv_files = [ legsta, mykrobe, ngmaster ]
tsv_files = [ drprg, checkm2, elgato, meningotype, seqsero2, seqsero2s, shigapass, kleborate, mlst, emmtyper, pbptyper ]

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

summary_df = pd.read_csv(names, dtype = str, names=input_cols, index_col=None, header=0, delimiter=",")
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

# fastqc : merging relevant rows into one
if exists(fastqc):
    file = fastqc
    print("Adding results for " + file)
    analysis = "fastqc"
    new_df = pd.read_csv(file, dtype=str, index_col=False)
    
    # FastQC output usually has two rows per sample (R1 and R2).
    # We split them and merge them side-by-side to get totals.
    R1_df = new_df.drop_duplicates(subset='sample', keep="first").add_prefix('R1_')
    R2_df = new_df.drop_duplicates(subset='sample', keep="last").add_prefix('R2_')
    
    new_df = pd.merge(R1_df, R2_df, left_on="R1_sample", right_on="R2_sample", how='left')
    new_df['sample'] = new_df['R1_sample']
    
    # Calculate totals
    new_df[['R1_Total Sequences', 'R2_Total Sequences', 'R1_Sequences flagged as poor quality', 'R2_Sequences flagged as poor quality']] = new_df[['R1_Total Sequences', 'R2_Total Sequences', 'R1_Sequences flagged as poor quality', 'R2_Sequences flagged as poor quality']].astype(int)
    new_df['total_sequences'] = new_df['R1_Total Sequences'] + new_df['R2_Total Sequences']
    new_df['flagged_sequences'] = new_df['R1_Sequences flagged as poor quality'] + new_df['R2_Sequences flagged as poor quality']
    new_df['percent_flagged'] = (new_df['flagged_sequences'] / new_df['total_sequences']) * 100

    # Rename and Merge
    new_df = new_df[['sample', 'total_sequences', 'flagged_sequences', 'percent_flagged']]
    new_df = new_df.add_prefix('fastqc_')
    
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on="fastqc_sample", how='left')
    summary_df.drop("fastqc_sample", axis=1, inplace=True)

# kaptive : merging relevant rows into one
if exists(kaptive) :
    file = kaptive
    print("Adding results for " + file)
    analysis = "kaptive"
    new_df = pd.read_table(file, dtype = str, index_col= False, sep="\t")
    new_df = new_df.add_prefix(analysis + '_')
    new_df.columns = [x.lower() for x in new_df.columns]
    K_df   = new_df[new_df['kaptive_best match locus'].str.contains("K")].copy()
    K_df   = K_df.add_suffix('_K')
    O_df   = new_df[new_df['kaptive_best match locus'].str.contains("O")].copy()
    O_df   = O_df.add_suffix('_O')

    summary_df = pd.merge(summary_df, O_df, left_on="sample", right_on=analysis + "_assembly_O", how = 'left')
    summary_df.drop(analysis + "_assembly_O", axis=1, inplace=True)
    summary_df = pd.merge(summary_df, K_df, left_on="sample", right_on=analysis + "_assembly_K", how = 'left')
    summary_df.drop(analysis + "_assembly_K", axis=1, inplace=True)

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
    new_df = new_df.replace('Sample', 'sample', regex=True)
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
if exists(serotypefinder) :
    file = serotypefinder
    print("Adding results for " + file)
    analysis = "serotypefinder"
    new_df = pd.read_table(file, dtype = str, index_col= False)
    new_df = new_df.sort_values(by='Identity', ascending=False)
    new_df = new_df.drop_duplicates(subset=['sample', 'Database'], keep="first")
    new_df = new_df.add_prefix(analysis + '_')
    H_df   = new_df[new_df[analysis + '_Database' ] == 'H_type'].copy()
    H_df   = H_df.add_suffix('_H')
    O_df   = new_df[new_df[analysis + '_Database' ] == 'O_type'].copy()
    O_df   = O_df.add_suffix('_O')
    summary_df = pd.merge(summary_df, O_df, left_on="sample", right_on=analysis + "_sample_O", how = 'left')
    summary_df.drop(analysis + "_sample_O", axis=1, inplace=True)
    summary_df = pd.merge(summary_df, H_df, left_on="sample", right_on=analysis + "_sample_H", how = 'left')
    summary_df.drop(analysis + "_sample_H", axis=1, inplace=True)


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
    if "FastQC_mqc-generalstats-fastqc-avg_sequence_length" in new_df.columns :
        tmp_df = new_df[["Sample","FastQC_mqc-generalstats-fastqc-avg_sequence_length"]].copy()
        tmp_df["fastqc_avg_length"] = tmp_df["FastQC_mqc-generalstats-fastqc-avg_sequence_length"]
        tmp_df.drop("FastQC_mqc-generalstats-fastqc-avg_sequence_length", axis=1, inplace=True)
        tmp_df = tmp_df.dropna(subset=['fastqc_avg_length'])
        
        summary_df["possible_fastqc_name"] = summary_df['file'].str.split(" ").str[0].str.split(".").str[0]
        summary_df = pd.merge(summary_df, tmp_df, left_on="possible_fastqc_name", right_on="Sample", how = 'left')
        summary_df.drop("Sample", axis=1, inplace=True)
        summary_df.drop("possible_fastqc_name", axis=1, inplace=True)

    if "fastp-pct_surviving" in new_df.columns :
        tmp_df = new_df[["Sample","fastp-pct_surviving"]].copy()
        tmp_df["fastp_pct_passed_reads"] = tmp_df["fastp-pct_surviving"].astype(float).round(2)
        tmp_df.drop("fastp-pct_surviving", axis=1, inplace=True)
        tmp_df = tmp_df.dropna(subset=['fastp_pct_passed_reads'])
        
        summary_df["possible_fastp_name"] = summary_df['file'].str.split(" ").str[0].str.split(".").str[0]
        summary_df = pd.merge(summary_df, tmp_df, left_on="possible_fastp_name", right_on="Sample", how = 'left')
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

# 1. Start with an empty column
summary_df['predicted_organism'] = pd.NA

# 2. Priority 1: Skani (Genome-level alignment)
if 'skani_organism' in summary_df.columns:
    summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['skani_organism'])

# 3. Priority 2: Kraken2 (Read-level k-mer analysis)
if 'kraken2_top_organism' in summary_df.columns:
    summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['kraken2_top_organism'])

# 4. Priority 3: Mash Dist (Fast k-mer sketching)
if 'mash_dist_organism' in summary_df.columns:
    summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['mash_dist_organism'])

# 5. Final Cleanup
summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna("Unknown")

# TODO : coverage estimates
# TODO : warnings
# TODO : E. coli / Shigella differentiation (add shigapass results)

# ##########################################
# # size and coverage estimates            #
# ##########################################

# if "fastqc_total_sequences" in summary_df.columns and 'fastqc_avg_length' in summary_df.columns:
#     print("Estimating coverage")

#     # 1. Calculate Total Raw Bases
#     summary_df['total_bases'] = summary_df['fastqc_total_sequences'].astype(float) * summary_df['fastqc_avg_length'].astype(float)
    
#     # 2. Reference-based Coverage (Priority 1)
#     if exists(genome_sizes):
#         # Load your JSON mapping: {"Escherichia coli": 5000000, ...}
#         with open(genome_sizes, 'r') as f:
#             size_dict = json.load(f)
        
#         # Map the predicted organism to its expected size
#         summary_df['expected_size'] = summary_df['predicted_organism'].map(size_dict)
#         summary_df['rep_estimated_coverage'] = summary_df['total_bases'] / summary_df['expected_size'].astype(float)

#     # 3. QUAST-based Coverage (Priority 2)
#     # Using 'total_length' because it's the actual size of your specific assembly
#     if 'quast_total_length' in summary_df.columns:
#         summary_df['quast_total_length'] = pd.to_numeric(summary_df['quast_total_length'], errors='coerce')
#         summary_df['quast_estimated_coverage'] = summary_df['total_bases'] / summary_df['quast_total_length']

#     # 4. Final Coverage Waterfall
#     # We start with Reference-based, then fill gaps with QUAST, then Mash
#     summary_df['final_coverage'] = pd.NA

#     if 'rep_estimated_coverage' in summary_df:
#         summary_df['final_coverage'] = summary_df['final_coverage'].fillna(summary_df['rep_estimated_coverage'])

#     if 'quast_estimated_coverage' in summary_df:
#         summary_df['final_coverage'] = summary_df['final_coverage'].fillna(summary_df['quast_estimated_coverage'])

#     if 'mash_err_mash_estimated_coverage' in summary_df:
#         summary_df['final_coverage'] = summary_df['final_coverage'].fillna(summary_df['mash_err_mash_estimated_coverage'].astype(float))

#     # Round for the final report
#     summary_df['final_coverage'] = pd.to_numeric(summary_df['final_coverage'], errors='coerce').round(2)



# if 'expected_size' in summary_df and 'quast_total_length' in summary_df:
#     summary_df['size_diff_ratio'] = (summary_df['quast_total_length'] / summary_df['expected_size']).astype(float)
    
#     # Flag if assembly is >20% different than expected
#     summary_df['size_warnings'] = summary_df['size_diff_ratio'].apply(
#         lambda x: "Assembly size differs from expected" if (x > 1.2 or x < 0.8) else ""
#     )
#     summary_df = add_warning(summary_df, 'size_warnings')

# # replacing Shigella with E. coli if ipaH+
# if 'predicted_organism' and 'shigatyper_hit' in summary_df.columns:
#     summary_df.loc[(summary_df['predicted_organism'].str.contains('Shigella')) & (~summary_df['shigatyper_hit'].str.contains('ipaH').notna()), 'predicted_organism'] = 'Escherichia coli'

# # genome size columns : checkm2_genome_size, quast_total_length_(>=_0_bp)



# # adding warnings for end user

#     amrfinder: for when big five are found

#     drprg?

#     el gato?

#     emmtyper?

#     fastp:

#     fastqc : 

#     kaptive?

#     kleborate:

#     kraken2:

#     menintotype?

#     mlst?

#     mykrobe?

#     ngmaster?

#     pbptyper?

#     plasmidfinder?

#     serotypefinder?

#     shigapass?



#     checkm2 : checkm2_completeness	checkm2_contamination

#     seqsero2 and seqsero2s: seqsero2_predicted_antigenic_profile	seqsero2_predicted_serotype

#     seqsero2 ST isn't the same as MLST ST

#     fastqc : fastqc_flagged_sequences	fastqc_percent_flagged


#     # Logic for warnings
#     def check_fastqc_warnings(row):
#         notes = []
#         if row['total_sequences'] < 10000: notes.append("Low read count (<10k)")
#         if row['percent_flagged'] > 1: notes.append("High % flagged sequences (>1%)")
#         return ", ".join(notes)

#     new_df['fastqc_warnings'] = new_df.apply(check_fastqc_warnings, axis=1)
    
#     # 5. Coverage Warnings
#     summary_df['coverage_warnings'] = summary_df['final_coverage'].apply(
#         lambda x: "Low coverage (<20x)" if x < 20 else ""
#     )
#     summary_df = add_warning(summary_df, 'coverage_warnings')


#         # --- 2. WARNING LOGIC ---
#     # Look for Carbapenemases specifically based on the Subclass column
#     carb_mask = amr_df['Subclass'].str.contains('CARBAPENEM', na=False, case=False)
#     amr_df['amr_warnings'] = ''
#     amr_df.loc[carb_mask, 'amr_warnings'] = amr_df.loc[carb_mask, 'Element symbol'].apply(
#         lambda x: f"Carbapenemase ({x})"
#     )

#     # Group warnings by Sample
#     amr_warn_summary = amr_df.groupby('Name')['amr_warnings'].apply(
#         lambda x: ", ".join(filter(None, x.unique()))
#     ).reset_index()
#     amr_warn_summary.columns = ['sample', 'amr_summary_warnings']

#     # Merge warnings into summary_df and use helper
#     summary_df = pd.merge(summary_df, amr_warn_summary, on='sample', how='left')
#     summary_df = add_warning(summary_df, 'amr_summary_warnings')
#     #print(summary_df)
# # to do, fix this
# # summary_df['warnings'] = summary_df['warnings'] + summary_df['kleborate_qc_warnings']

#     # adding warning logic

#     counts['warnings'] = counts['predicted_organisms'].apply(
#         lambda x: "mash dist predicted >10 organisms" if x > 10 else ""
#     )

#         counts['warnings'] = counts['predicted_organisms'].apply(
#         lambda x: "mash screen predicted >10 organisms" if x > 10 else ""
#     )
#     # Warning if more than 2 organisms are found at significant levels
#     grouped['kraken2_warnings'] = grouped['report'].apply(lambda x: "Kraken2 detects multiple species" if len(x) > 2 else "")

#     new_df['quast_warnings'] = ''

#     checks = {
#         'quast_n50': (25000, 'low N50', '<'),
#         'quast_l50': (100, 'high L50', '>'),
#         'quast_avg. coverage depth': (30, 'low coverage', '<'),
#         "quast_# n's per 100 kbp": (500, 'high gap content (>500 Ns/100kbp)', '>')
#     }

#     for col, (thresh, msg, op) in checks.items():
#         if col in new_df.columns:
#             new_df[col] = pd.to_numeric(new_df[col], errors='coerce')
#             mask = (new_df[col] < thresh) if op == '<' else (new_df[col] > thresh)
#             new_df.loc[mask, 'quast_warnings'] = new_df.loc[mask, 'quast_warnings'].apply(
#                 lambda x: f"{x}, {msg}" if x else msg
#             )
#     counts['skani_warnings'] = counts['organism'].apply(lambda x: "Skani predicted >10 organisms" if x > 10 else "")

#     # Sort and create report
#     purity_warnings = new_df.sort_values(['sample', 'Taxonomic_abundance'], ascending=[True, False])
    
#     # Now this line will work because 'check_isolate_purity' is defined at the top!
#     purity_report = purity_warnings.groupby('sample').apply(check_isolate_purity).reset_index()
#     purity_report.columns = ['sample', 'sylph_warnings']

#     summary_df['core_genome_warnings'] = summary_df['per_core_genome_genes'].apply(lambda x: "Low core genes," if x <= 85 else "")
#     summary_df['warnings']            = summary_df['warnings'] + summary_df['core_genome_warnings']


##########################################
# creating files                         #
##########################################

summary_df = summary_df.sort_values(by='sample')
summary_df = summary_df.fillna("")

summary_df.columns = summary_df.columns.str.replace(' ', '_')

summary_df.to_csv(extended + '.tsv', index=False, sep="\t")
summary_df.to_csv(extended + '.txt', index=False, sep=";")

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
    'quast_gc_(%)',
    'warnings',
    'amrfinder_genes_(per_cov/per_ident)',

    # species
    'predicted_organism',
    'mlst_matching_pubmlst_scheme',
    'mlst_st',
    'fastani_top_organism',
    'fastani_top_reference',
    'fastani_top_ani_estimate',
    'fastani_top_total_query_sequence_fragments',
    'fastani_top_fragments_aligned_as_orthologous_matches',
    'mash_reference',
    'mash_mash-distance',
    'mash_p-value',
    'mash_matching-hashes',
    'mash_organism',
    'plasmidfinder_plasmid_(identity)',

    # contamination
    # add sylph
    'kraken2_organism_(per_fragment)',

    # species specific information
    'seqsero2_predicted_antigenic_profile',
    'seqsero2_predicted_serotype',
    'emmtyper_predicted_emm-type',
    'kleborate_virulence_score',
    'kleborate_resistance_score',
    'kaptive_best_match_locus_O',
    'kaptive_best_match_locus_K',
    'elgato_st',
    'meningotype_serogroup',
    'mykrobe_phylo_group',
    'mykrobe_species',
    'mykrobe_lineage',
    'drprg_susceptibility',
    'pbptyper_pbptype',
    'serotypefinder_Serotype_O',
    'serotypefinder_Serotype_H'
    ]

set_columns = []
for new_column in final_columns :
    if new_column in summary_df.columns :
        set_columns.append(new_column)

summary_df.to_csv(final + '.tsv', columns = ['sample','file','version'] + set_columns, index=False, sep="\t")
summary_df.to_csv(final + '.txt', columns = ['sample','file','version'] + set_columns, index=False, sep=";")
