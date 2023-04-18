#!/bin/python

##########################################
# written by Erin Young                  #
# for creating summary file for grandeur #
##########################################

import pandas as pd
from os.path import exists

##########################################
# defining files                         #
##########################################

# input files
names          = 'input_files.csv'
amrfinderplus  = 'amrfinderplus.txt'
blobtools      = 'blobtools_summary.txt'
datasets       = 'datasets_summary.csv'
emmtyper       = 'emmtyper_summary.tsv'
fastani        = 'fastani_summary.csv'
fastqc         = 'fastqc_summary.csv'
fastqscan      = 'fastqscan_summary.csv'
kleborate      = 'kleborate_results.tsv'
kraken2        = 'kraken2_summary.csv'
legsta         = 'legsta_summary.csv'
mash           = 'mash_summary.csv'
mlst           = 'mlst_summary.tsv'
pbptyper       = 'pbptyper_summary.tsv'
plasmidfinder  = 'plasmidfinder_result.tsv'
quast          = 'quast_report.tsv'
seqsero2       = 'seqsero2_results.txt'
serotypefinder = 'serotypefinder_results.txt'
shigatyper     = 'shigatyper_results.txt'

# final files
final          = 'grandeur_summary'
extended       = 'grandeur_extended_summary'

##########################################
# grouping similar files                 #
##########################################

csv_files = [ fastqscan, legsta ]
tsv_files = [ quast, seqsero2, kleborate, mlst, emmtyper , pbptyper]

top_hit    = [ fastani ]

##########################################
# creating the summary dataframe         #
##########################################

summary_df = pd.read_csv(names, dtype = str)
columns = list(summary_df.columns)

# csv files
for file in csv_files :
    if exists(file) :
        print("Adding results for " + file)
        analysis = str(file).split("_")[0]
        new_df = pd.read_csv(file, dtype = str, index_col= False)
        new_df = new_df.add_prefix(analysis + "_")
        summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
        summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# csv files
for file in tsv_files :
    if exists(file) :
        print("Adding results for " + file)
        analysis = str(file).split("_")[0]
        new_df = pd.read_table(file, dtype = str, index_col= False)
        new_df = new_df.add_prefix(analysis + "_")
        summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
        summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# extracting top hit for multiline result
for file in top_hit :
    if exists(file) : 
        print("Adding results for " + file)
        analysis = str(file).split("_")[0]
        new_df = pd.read_csv(file, dtype = str, index_col= False)
        new_df = new_df.drop_duplicates(subset=['sample'], keep='first')
        new_df = new_df.add_prefix(analysis + "_")
        summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
        summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# for specific tools

# amrfinderplus : merging many rows into one with relevant information
if exists(amrfinderplus) :
    file = amrfinderplus
    print("Adding results for " + file)
    analysis = "amrfinder"
    new_df = pd.read_table(file, dtype = str, index_col= False)
    new_df['genes (per cov/per ident)'] = new_df['Gene symbol'] + ' (' + new_df['% Coverage of reference sequence'] + '/' + new_df['% Identity to reference sequence'] + ')'
    new_df = new_df[['Name', 'genes (per cov/per ident)']]
    new_df = new_df.groupby('Name', as_index=False).agg({'genes (per cov/per ident)': lambda x: list(x)})
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_Name", how = 'left')
    summary_df.drop(analysis + "_Name", axis=1, inplace=True)

# blobtools : merging many rows into one with relevant information
if exists(blobtools) :
    file = blobtools
    print("Adding results for " + file)
    analysis = "blobtools"
    new_df = pd.read_table(file, dtype = str, index_col= False)
    new_df = new_df[new_df['name'] != 'all']
    new_df['organism (per mapped reads)'] = new_df['name'] + ' (' + new_df['bam0_read_map_p'] + ')'
    new_df = new_df[['sample', 'organism (per mapped reads)']]
    new_df = new_df.groupby('sample', as_index=False).agg({'organism (per mapped reads)': lambda x: list(x)})
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# fastani : merging relevant rows into one
if exists(fastani) :
    file = fastani
    print("Adding results for " + file)
    analysis = "fastani"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    new_df['genome (ANI estimate)'] = new_df['reference'].str.split('_').str[0] + " " + new_df['reference'].str.split("_").str[1] + " " + new_df['reference'].str.split('_').str[-2] + "_" + new_df['reference'].str.split('_').str[-1] + " (" + new_df['ANI estimate'] + ")"
    new_df = new_df[['sample', 'genome (ANI estimate)']]
    new_df = new_df.groupby('sample', as_index=False).agg({'genome (ANI estimate)': lambda x: list(x)})
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# fastqc : merging relevant rows into one
if exists(fastqc) :
    file = fastqc
    print("Adding results for " + file)
    analysis = "fastqc"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    R1_df = new_df.drop_duplicates(subset='sample', keep="first")
    R1_df = R1_df.add_prefix('R1_')
    R2_df = new_df.drop_duplicates(subset='sample', keep="last")
    R2_df = R2_df.add_prefix('R2_')
    new_df = pd.merge(R1_df, R2_df, left_on="R1_sample", right_on="R2_sample", how = 'left')
    new_df['sample'] = new_df['R1_sample']
    new_df[['R1_Total Sequences', 'R2_Total Sequences', 'R1_Sequences flagged as poor quality', 'R2_Sequences flagged as poor quality']] = new_df[['R1_Total Sequences', 'R2_Total Sequences', 'R1_Sequences flagged as poor quality', 'R2_Sequences flagged as poor quality']].astype(int)
    new_df['total sequences']    = new_df.apply(lambda x: x['R1_Total Sequences'] + x['R2_Total Sequences'], axis=1)
    new_df['flagged sequences']  = new_df.apply(lambda x: x['R1_Sequences flagged as poor quality'] + x['R2_Sequences flagged as poor quality'], axis=1) 
    new_df = new_df[['sample','total sequences', 'flagged sequences']]
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# kraken2 : merging relevant rows into one
if exists(kraken2) :
    file = kraken2
    print("Adding results for " + file)
    analysis = "kraken2"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    new_df['organism (per fragment)'] = new_df['Scientific name'] + " (" + new_df['Percentage of fragments'] + ' ' + new_df['Type'] + ")"    
    new_df = new_df[['Sample', 'organism (per fragment)']]
    new_df = new_df.groupby('Sample', as_index=False).agg({'organism (per fragment)': lambda x: list(x)})
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_Sample", how = 'left')
    summary_df.drop(analysis + "_Sample", axis=1, inplace=True)

# mash : top hit of potentially two different files
if exists(mash) :
    file = mash
    print("Adding results for " + file)
    analysis = "mash"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    new_df.sort_values(by = ['P-value', 'mash-distance'], ascending = [True, True], inplace=True)
    new_df = new_df.drop_duplicates(subset=['sample'], keep='first')
    new_df = new_df.add_prefix(analysis + "_")
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# plasmidfinder : merging relevant rows into one
if exists(plasmidfinder) :
    file = plasmidfinder
    print("Adding results for " + file)
    analysis = "plasmidfinder"
    new_df = pd.read_table(file, dtype = str, index_col= False)
    new_df['plasmid (identity)'] = new_df['Plasmid'] + " (" + new_df['Identity'] + ")"
    new_df = new_df[['sample', 'plasmid (identity)']]
    new_df = new_df.groupby('sample', as_index=False).agg({'plasmid (identity)': lambda x: list(x)})
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
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
    H_df   = new_df[new_df[analysis + '_Database' ] == 'H_type']
    H_df   = H_df.add_suffix('_H')
    O_df   = new_df[new_df[analysis + '_Database' ] == 'O_type']
    O_df   = O_df.add_suffix('_O')
    summary_df = pd.merge(summary_df, O_df, left_on="sample", right_on=analysis + "_sample_O", how = 'left')
    summary_df.drop(analysis + "_sample_O", axis=1, inplace=True)
    summary_df = pd.merge(summary_df, H_df, left_on="sample", right_on=analysis + "_sample_H", how = 'left')
    summary_df.drop(analysis + "_sample_H", axis=1, inplace=True)

# shigatyper : combining rows into one
if exists(shigatyper) :
    file = shigatyper
    print("Adding results for " + file)
    analysis = "shigatyper"
    new_df = pd.read_table(file, dtype = str, index_col= False)
    new_df = new_df.groupby('sample', as_index=False).agg({'Hit': lambda x: list(x)})
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

##########################################
# creating files                         #
##########################################

summary_df.columns = summary_df.columns.str.replace(' ', '_')

summary_df.to_csv(extended + '.tsv', index=False, sep="\t")
summary_df.to_csv(extended + '.txt', index=False, sep=";")

# reducing to the top 1 or 2 results for each analysis
final_columns = [
# general information
'fastqc_total_sequences',
'fastqc_flagged_sequences',
'fastqscan_coverage',
'quast_#_contigs',
'quast_GC_(%)',
'amrfinder_genes_(per_cov/per_ident)',

# species
'mlst_matching_PubMLST_scheme',
'mlst_ST',
'mash_reference',
'mash_mash-distance',
'mash_P-value',
'mash_matching-hashes',
'mash_organism',
'fastani_reference',
'fastani_ANI_estimate',
'fastani_total_query_sequence_fragments',
'fastani_fragments_aligned_as_orthologous_matches',
'plasmidfinder_plasmid_(identity)',

# contamination
'blobtools_organism_(per_mapped_reads)',
'kraken2_organism_(per_fragment)',

# species specific information
'seqsero2_Predicted_antigenic_profile',
'seqsero2_Predicted_serotype',
'emmtyper_Predicted_emm-type',
'kleborate_virulence_score',
'pbptyper_pbptype',
'legsta_SBT',
'serotypefinder_Serotype_O',
'serotypefinder_Serotype_H',
'shigatyper_Hit']

set_columns = []
for new_column in final_columns :
    if new_column in summary_df.columns :
        set_columns.append(new_column)

summary_df.to_csv(final + '.tsv', columns = ['sample','file','version','reads','phix_reads'] + set_columns, index=False, sep="\t")
summary_df.to_csv(final + '.txt', columns = ['sample','file','version','reads','phix_reads'] + set_columns, index=False, sep=";")