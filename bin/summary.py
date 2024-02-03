#!/bin/python

##########################################
# written by Erin Young                  #
# for creating summary file for grandeur #
##########################################

import pandas as pd
import json
from os.path import exists

##########################################
# defining files                         #
##########################################

# input files
names          = 'input_files.csv'
amrfinderplus  = 'amrfinderplus.txt'
blobtools      = 'blobtools_summary.txt'
core           = 'multiqc_core_genome_evaluation-plot.txt'
datasets       = 'datasets_summary.csv'
drprg          = 'replaceme.csv'
elgato         = 'replaceme.csv'
emmtyper       = 'emmtyper_summary.tsv'
fastani        = 'fastani_summary.csv'
fastani_len    = "fastani_top_len.csv"
fastqc         = 'fastqc_summary.csv'
genome_sizes   = "genome_sizes.json"
kleborate      = 'kleborate_results.tsv'
kraken2        = 'kraken2_summary.csv'
legsta         = 'legsta_summary.csv'
mash           = 'mash_summary.csv'
mash_err       = 'mash_err_summary.csv'
mlst           = 'mlst_summary.tsv'
mykrobe        = 'mykrobe_summary.csv'
pbptyper       = 'pbptyper_summary.tsv'
plasmidfinder  = 'plasmidfinder_result.tsv'
quast          = 'quast_report.tsv'
seqsero2       = 'seqsero2_results.txt'
serotypefinder = 'serotypefinder_results.txt'
shigatyper     = 'shigatyper_results.txt'
multiqc_json   = 'multiqc_data.json'
multiqc_stats  = 'multiqc_general_stats.txt'

# final files
final          = 'grandeur_summary'
extended       = 'summary/grandeur_extended_summary'

##########################################
# grouping similar files                 #
##########################################

csv_files = [ legsta, mykrobe ]
tsv_files = [ quast, seqsero2, kleborate, mlst, emmtyper , pbptyper]
top_hit   = [ fastani ]

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

summary_df = pd.read_csv(names, dtype = str, index_col=None)
summary_df['warnings'] = ''
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
    new_df = new_df.sort_values('Gene symbol')
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
    new_df = new_df.sort_values('bam0_read_map_p', ascending = False)
    
    tmp_df = new_df.drop_duplicates(subset='sample', keep="first").copy()
    tmp_df = tmp_df[['sample', 'name']]
    tmp_df['top_organism'] = tmp_df['name']
    tmp_df = tmp_df.add_prefix(analysis + '_')
    
    new_df['organism (per mapped reads)'] = new_df['name'] + ' (' + new_df['bam0_read_map_p'] + ')'
    new_df = new_df[['sample', 'organism (per mapped reads)']]
    new_df = new_df.groupby('sample', as_index=False).agg({'organism (per mapped reads)': lambda x: list(x)})
    new_df['warning'] = new_df['organism (per mapped reads)'].apply(lambda x: "Blobtools multiple organisms," if ','.join(x).count(',') >= 2 else "")
    new_df = new_df.add_prefix(analysis + '_')
    
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)
    summary_df = pd.merge(summary_df, tmp_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df['warnings'] = summary_df['warnings'] + summary_df['blobtools_warning']

# fastani : merging relevant rows into one
if exists(fastani) :
    file = fastani
    print("Adding results for " + file)
    analysis = "fastani"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    new_df['genome (ANI estimate)'] = new_df['reference'].str.split('_').str[0] + " " + new_df['reference'].str.split("_").str[1] + " " + new_df['reference'].str.split('_').str[-2] + "_" + new_df['reference'].str.split('_').str[-1] + " (" + new_df['ANI estimate'] + ")"
    new_df = new_df[['sample', 'genome (ANI estimate)']]
    new_df = new_df.groupby('sample', as_index=False).agg({'genome (ANI estimate)': lambda x: list(x)})
    new_df['warning'] = new_df['genome (ANI estimate)'].apply(lambda x: "Multiple FastANI hits," if ','.join(x).count(',') >= 4 else "")
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)
    summary_df['warnings'] = summary_df['warnings'] + summary_df['fastani_warning']

# fastqc : merging relevant rows into one
if exists(fastqc) :
    file = fastqc
    print("Adding results for " + file)
    analysis = "fastqc"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    R1_df = new_df.drop_duplicates(subset='sample', keep="first").copy()
    R1_df = R1_df.add_prefix('R1_')
    R2_df = new_df.drop_duplicates(subset='sample', keep="last").copy()
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
    new_df = new_df.sort_values(['Sample', 'Percentage of fragments'], ascending= False)
    new_df['top_organism'] = new_df['Scientific name']

    tmp_df = new_df.drop_duplicates(subset=['Sample'], keep="first").copy()
    tmp_df = tmp_df[['Sample', 'top_organism']]
    tmp_df = tmp_df.add_prefix(analysis + '_')

    new_df['organism (per fragment)'] = new_df['Scientific name'] + " (" + new_df['Percentage of fragments'] + ")"    
    new_df = new_df[['Sample', 'organism (per fragment)']]
    new_df = new_df.groupby('Sample', as_index=False).agg({'organism (per fragment)': lambda x: list(x)})
    new_df['warning'] = new_df['organism (per fragment)'].apply(lambda x: "Kraken2 multiple organisms," if ','.join(x).count(',') >= 2 else "")
    new_df = new_df.add_prefix(analysis + '_')
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_Sample", how = 'left')
    summary_df.drop(analysis + "_Sample", axis=1, inplace=True)
    summary_df = pd.merge(summary_df, tmp_df, left_on="sample", right_on=analysis + "_Sample", how = 'left')
    summary_df.drop(analysis + "_Sample", axis=1, inplace=True)
    summary_df['warnings'] = summary_df['warnings'] + summary_df['kraken2_warning']

# mash : top hit of potentially two different files
if exists(mash) :
    file = mash
    print("Adding results for " + file)
    analysis = "mash"
    new_df = pd.read_csv(file, dtype = str, index_col= False)
    new_df = new_df.sort_values(by = ['P-value', 'mash-distance'], ascending = [True, True])
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
    H_df   = new_df[new_df[analysis + '_Database' ] == 'H_type'].copy()
    H_df   = H_df.add_suffix('_H')
    O_df   = new_df[new_df[analysis + '_Database' ] == 'O_type'].copy()
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

# multiqc : bbduk and fastp
if exists(multiqc_json) :
    file = multiqc_json
    print("Adding analysis parsed via multiqc in " + file)
    with open(file) as multiqc_data:
        data = json.load(multiqc_data)
        
        # fastp filtered reads
        if "fastp_filtered_reads_plot" in data["report_plot_data"].keys():
            samples = [sample.replace("_rmphix_R1", "") for sample in data["report_plot_data"]['fastp_filtered_reads_plot']['samples'][0]]
            fastp_passedreads_df = pd.DataFrame(samples, columns=['fastp_sample'])
            fastp_passedreads_df['fastp_passed_reads'] = data["report_plot_data"]['fastp_filtered_reads_plot']['datasets'][0][0]['data']
            summary_df = pd.merge(summary_df, fastp_passedreads_df, left_on="sample", right_on="fastp_sample", how = 'left')
            summary_df.drop("fastp_sample", axis=1, inplace=True)
        
        # bbduk phix reads
        if "bbmap" in data['report_saved_raw_data'].keys():
            print("Adding in phix reads from bbmap")
            samples = [sample.replace(".phix", "") for sample in data['report_saved_raw_data']['bbmap']['stats'].keys()]
            phix_reads=[]
            for sample in data['report_saved_raw_data']['bbmap']['stats'].keys() :
                phix_reads.append(data['report_saved_raw_data']['bbmap']['stats'][sample]['kv']['Matched'])
            bbduk_phixreads_df = pd.DataFrame(samples, columns=['bbduk_sample'])
            bbduk_phixreads_df['bbduk_phix_reads'] = phix_reads
            summary_df = pd.merge(summary_df, bbduk_phixreads_df, left_on="sample", right_on="bbduk_sample", how = 'left')
            summary_df.drop("bbduk_sample", axis=1, inplace=True)

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

    summary_df['core_genome_warning'] = summary_df['per_core_genome_genes'].apply(lambda x: "Low core genes," if x <= 85 else "")
    summary_df['warnings']            = summary_df['warnings'] + summary_df['core_genome_warning']


# size : getting the size and coverage and warning if there's too much stdev
        
# getting different size estimates
# from circulocov

##########################################
# predicting organism                    #
##########################################

print("Predicting organism")

summary_df['predicted_organism'] = pd.NA

# TODO : add in fastani

if 'blobtools_top_organism' in summary_df:
    def fill_predicted_organism(row):
        if pd.isna(row['predicted_organism']):
            return row['blobtools_top_organism']
        else:
            return row['predicted_organism']

    summary_df['predicted_organism'] = summary_df.apply(fill_predicted_organism, axis=1)

if 'kraken2_top_organism' in summary_df:
    def fill_predicted_organism(row):
        if pd.isna(row['predicted_organism']):
            return row['kraken2_top_organism']
        else:
            return row['predicted_organism']

    summary_df['predicted_organism'] = summary_df.apply(fill_predicted_organism, axis=1)

if 'mash_organism' in summary_df:
    def fill_predicted_organism(row):
        if pd.isna(row['predicted_organism']):
            return row['mash_organism']
        else:
            return row['predicted_organism']

    summary_df['predicted_organism'] = summary_df.apply(fill_predicted_organism, axis=1)

##########################################
# size and coverage estimates            #
##########################################

print("Estimating coverage")

if "fastqc_total sequences" and 'fastqc_avg_length' in summary_df:
    print("Estimating coverage")

    summary_df['total_bases'] = summary_df['fastqc_total sequences'].astype(float) * summary_df['fastqc_avg_length'].astype(float)
    summary_df['coverage_for_1.5M_genome'] = summary_df['total_bases'].astype(float) /  1500000
    summary_df['coverage_for_2M_genome']   = summary_df['total_bases'].astype(float) /  2000000
    summary_df['coverage_for_2.5M_genome'] = summary_df['total_bases'].astype(float) /  2500000
    summary_df['coverage_for_3M_genome']   = summary_df['total_bases'].astype(float) /  3000000
    summary_df['coverage_for_3.5M_genome'] = summary_df['total_bases'].astype(float) /  3500000
    summary_df['coverage_for_4M_genome']   = summary_df['total_bases'].astype(float) /  4000000
    summary_df['coverage_for_4.5M_genome'] = summary_df['total_bases'].astype(float) /  4500000
    summary_df['coverage_for_5M_genome']   = summary_df['total_bases'].astype(float) /  5000000
    summary_df['coverage_for_5.5M_genome'] = summary_df['total_bases'].astype(float) /  5500000
    summary_df['coverage_for_6M_genome']   = summary_df['total_bases'].astype(float) /  6000000
    summary_df['coverage_for_6.5M_genome'] = summary_df['total_bases'].astype(float) /  6500000
    summary_df['coverage_for_7M_genome']   = summary_df['total_bases'].astype(float) /  7000000
    summary_df['coverage_for_7.5M_genome'] = summary_df['total_bases'].astype(float) /  7500000
    summary_df['coverage_for_8M_genome']   = summary_df['total_bases'].astype(float) /  8000000
    summary_df['coverage_for_8.5M_genome'] = summary_df['total_bases'].astype(float) /  8500000
    summary_df['coverage_for_9M_genome']   = summary_df['total_bases'].astype(float) /  9000000
    summary_df['coverage_for_9.5M_genome'] = summary_df['total_bases'].astype(float) /  9500000
    summary_df['coverage_for_10M_genome']  = summary_df['total_bases'].astype(float) / 10000000

    if exists(mash_err) :
        file = mash_err
        print("Adding cov and size estimates from " + file)
        new_df = pd.read_csv(file, dtype = str, index_col= False)
        new_df["mash_estimated_genome_size"] = new_df["mash_estimated_genome_size"].astype(float)
        new_df["mash_estimated_coverage"] = new_df["mash_estimated_coverage"].astype(float)
        summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on="sample", how = 'left')

    if exists(fastani_len):
        file = fastani_len
        print("getting size estimate from fastani top hit" + file)
        new_df = pd.read_csv(file, dtype = str, index_col= False)
        new_df['predicted_organism'] = new_df["top_hit"].astype(str).apply(lambda x: "_".join(x.split('_')[:2]))
        new_df.rename(columns={'top_len': 'fastani_estimated_genome_size'}, inplace=True) 
        summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on="sample", how = 'left')
        summary_df['fastani_estimated_coverage'] = summary_df['total_bases'].astype(float) / summary_df['fastani_estimated_genome_size'].astype(float)

    if exists(genome_sizes) :
        genome_size_df = pd.read_json(genome_sizes, orient='records')
        genome_size_df = genome_size_df.reset_index()
        genome_size_df.rename(columns={'genome_sizes': 'rep_estimated_genome_size', 'index': 'predicted_organism'}, inplace=True)

        summary_df = pd.merge(summary_df, genome_size_df, left_on='predicted_organism', right_on='predicted_organism', how = 'left')
        summary_df['rep_estimated_coverage'] = summary_df['total_bases'].astype(float) / summary_df['rep_estimated_genome_size'].astype(float)

    if 'quast_Total length' in summary_df:
        summary_df['quast_estimated_genome_size'] = summary_df['quast_Total length']
        summary_df['quast_estimated_coverage']    = summary_df['total_bases'].astype(float) / summary_df['quast_estimated_genome_size'].astype(float)

    summary_df['coverage'] = pd.NA

    cov_columns = []

    if 'rep_estimated_coverage' in summary_df:
        cov_columns.append('rep_estimated_coverage')
        def fill_coverage(row):
            if pd.isna(row['coverage']):
                return row['rep_estimated_coverage']
            else:
                return row['coverage']

        summary_df['coverage'] = summary_df.apply(fill_coverage, axis=1)

    if 'fastani_estimated_coverage' in summary_df:
        cov_columns.append('fastani_estimated_coverage')
        def fill_coverage(row):
            if pd.isna(row['coverage']):
                return row['fastani_estimated_coverage']
            else:
                return row['coverage']

        summary_df['coverage'] = summary_df.apply(fill_coverage, axis=1)

    if 'quast_estimated_coverage' in summary_df:
        cov_columns.append('quast_estimated_coverage')
        def fill_coverage(row):
            if pd.isna(row['coverage']):
                return row['quast_estimated_coverage']
            else:
                return row['coverage']

        summary_df['coverage'] = summary_df.apply(fill_coverage, axis=1)

    if 'mash_estimated_coverage' in summary_df:
        cov_columns.append('mash_estimated_coverage')
        def fill_coverage(row):
            if pd.isna(row['coverage']):
                return row['mash_estimated_coverage']
            else:
                return row['coverage']

        summary_df['coverage'] = summary_df.apply(fill_coverage, axis=1)

    if 'coverage_for_5M_genome' in summary_df:
        def fill_coverage(row):
            if pd.isna(row['coverage']):
                return row['coverage_for_5M_genome']
            else:
                return row['coverage']

        summary_df['coverage'] = summary_df.apply(fill_coverage, axis=1)

    summary_df['coverage'] = summary_df['coverage'].round(2)

    summary_df['coverage_warning'] = summary_df['coverage'].apply(lambda x: "Low coverage," if x <= 20 else "")
    summary_df['warnings']         = summary_df['warnings'] + summary_df['coverage_warning']

    if cov_columns:
        summary_df['cov_average'] = summary_df[cov_columns].mean(axis = 1, skipna = True)
        summary_df['cov_stdev'] = summary_df[cov_columns].std(axis = 1, skipna = True)
        summary_df['cov_stdev_ratio'] = summary_df['cov_stdev']/summary_df['cov_average']
        summary_df['warning'] = summary_df['cov_stdev_ratio'].apply(lambda x: "Variable genome size," if x >= 0.3 else "")

##########################################
# creating files                         #
##########################################

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
'fastp_passed_reads',
'bbduk_phix_reads',
'quast_#_contigs',
'quast_GC_(%)',
'warnings',
'amrfinder_genes_(per_cov/per_ident)',

# species
'predicted_organism',
'mlst_matching_PubMLST_scheme',
'mlst_ST',
'fastani_reference',
'fastani_ANI_estimate',
'fastani_total_query_sequence_fragments',
'fastani_fragments_aligned_as_orthologous_matches',
'mash_reference',
'mash_mash-distance',
'mash_P-value',
'mash_matching-hashes',
'mash_organism',
'plasmidfinder_plasmid_(identity)',

# contamination
'blobtools_organism_(per_mapped_reads)',
'kraken2_organism_(per_fragment)',

# species specific information
'seqsero2_Predicted_antigenic_profile',
'seqsero2_Predicted_serotype',
'emmtyper_Predicted_emm-type',
'kleborate_virulence_score',
'kleborate_resistance_score',
'mykrobe_phylo_group',
'mykrobe_species',
'mykrobe_lineage',
'pbptyper_pbptype',
'legsta_SBT',
'serotypefinder_Serotype_O',
'serotypefinder_Serotype_H',
'shigatyper_Hit']

set_columns = []
for new_column in final_columns :
    if new_column in summary_df.columns :
        set_columns.append(new_column)

summary_df.to_csv(final + '.tsv', columns = ['sample','file','version'] + set_columns, index=False, sep="\t")
summary_df.to_csv(final + '.txt', columns = ['sample','file','version'] + set_columns, index=False, sep=";")