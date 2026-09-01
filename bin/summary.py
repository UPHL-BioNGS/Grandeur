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

from bin.summary_parsers import parse_concatenated_json, parse_file
from bin.summary_phylo_utils import get_snp_distance_stats, get_tip_distance_stats
from bin.summary_qc_warnings import (
    add_checkm2_warnings,
    add_fastqc_warnings,
    add_quast_warnings,
    add_st_mismatch_warning,
    add_sylph_warnings,
    add_taxonomy_warnings,
)
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
from bin.summary_phylogenetics import summarize_phylogenetics
from bin.summary_plasmidfinder import summarize_plasmidfinder
from bin.summary_qc_warnings import add_qc_warnings
from bin.summary_quast import summarize_quast, summarize_quast_contig
from bin.summary_serotypefinder import summarize_serotypefinder
from bin.summary_skani import summarize_skani
from bin.summary_spestimator import summarize_spestimator
from bin.summary_sylph import summarize_sylph
from bin.summary_warnings import add_warnings

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

# establish coverage

# predict organism

# adding warnings

# create final files






































