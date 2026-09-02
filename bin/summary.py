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

from bin.summary_parsers import parse_file
from bin.summary_phylo_utils import get_snp_distance_stats, get_tip_distance_stats
from bin.summary_amrfinder import summarize_amrfinder
from bin.summary_core_genome import summarize_core_genome
from bin.summary_coverage import summarize_coverage
from bin.summary_create_files import create_final_summary, create_extended_summary
from bin.summary_datasets import summarize_datasets
from bin.summary_fastqc import summarize_fastqc
from bin.summary_kraken2 import summarize_kraken2
from bin.summary_mash_err import summarize_mash_err
from bin.summary_mash import summarize_mashdist, summarize_mashscreen
from bin.summary_multiqc import summarize_multiqc
from bin.summary_organism import predict_organism
from bin.summary_phylogenetics import summarize_newick, summarize_snpdist, summarize_gotree
from bin.summary_plasmidfinder import summarize_plasmidfinder
from bin.summary_quast import summarize_quast, summarize_quast_contig
from bin.summary_serotypefinder import summarize_serotypefinder
from bin.summary_skani import summarize_skani
from bin.summary_spestimator import summarize_spestimator
from bin.summary_sylph import summarize_sylph
from bin.summary_create_files import create_extended_summary, create_final_summary

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

#phylogenetics
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

if exists(amrfinderplus):
    summary_df = summarize_amrfinder(summary_df, amrfinderplus)

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
    summary_df = summarize_amrfinder(summary_df, amrfinderplus)

# datasets : adding a count of reference genomes available
if exists(datasets):
    summary_df = summarize_datasets(summary_df, datasets)

# fastqc
if exists(fastqc):
    summary_df = summarize_fastqc(summary_df, fastqc)

# kraken2 : merging relevant rows into one
if exists(kraken2):
    summary_df = summarize_kraken2(summary_df, kraken2)

# mash dist
if exists(mash_dist) :
    summary_df = summarize_mashdist(summary_df, mash_dist)

# mash screen
if exists(mash_screen) :
    summary_df = summarize_mashscreen(summary_df, mash_screen)

# getting genome size from mash err
if exists(mash_err):
    summary_df = summarize_mash_err(summary_df, mash_err)

# plasmidfinder : merging relevant rows into one
if exists(plasmidfinder) :
    summary_df = summarize_plasmidfinder(summary_df, plasmidfinder)

# quast : combining both files
q_df  = pd.DataFrame()
qc_df = pd.DataFrame()
if exists(quast):
    q_df = summarize_quast_reads(summary_df, quast)

if exists(quast_contig):
    qc_df = summarize_quast_contig(summary_df, quast_contig)

if exists(quast) or exists(quast_contig):
    summary_df = summarize_quast(summary_df, q_df, qc_df)

# serotypefinder : splitting O and H groups, getting the top hit for O and H group, combining rows
if exists(serotypefinder):
    summary_df = summarize_serotypefinder(summary_df, serotypefinder)

# skani
if exists(skani):
    summary_df = summarize_skani(summary_df, skani)

# spestimator : counting unique reference hits per sample
if exists(spestimator):
    summary_df = summarize_spestimator(summary_df, spestimator)

# sylph
if exists(sylph):
    summary_df = summarize_sylph(summary_df, sylph)

if exists(multiqc_stats) : 
    summary_df = summarize_multiqc(summary_df, multiqc_stats)

# core genome analysis file is also from multiqc
if exists(core):
    summary_df = summarize_core_genome(summary_df, core)

##########################################
# predicting organism                    #
##########################################

summary_df = predict_organism(summary_df)

##########################################
# size and coverage estimates            #
##########################################

summary_df = summarize_coverage(summary_df, genome_sizes)

##########################################
# summarizing phylogenetics              #
##########################################

for nwk in newick_files:
    if exists(nwk) :
        summary_df = summarize_phylogenetics(summary_df, nwk)

for snp_matrix in snpdist_matrices:
    if exists(snp_matrix) :
        summary_df = summarize_phylogenetics(summary_df, snp_matrix)

if exists(gotree):
    summary_df = summarize_phylogenetics(summary_df, gotree)



##########################################
# adding final flags and warnings        #
##########################################




##########################################
# creating files                         #
##########################################

summary_df = create_extended_summary(summary_df)
create_final_summary(summary_df)