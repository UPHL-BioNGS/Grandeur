#!/usr/bin/env python

import shutil
from os.path import exists
import pandas as pd
import numpy as np

##########################################
# defining files                         #
##########################################

# input files
amrfinder_input         = 'amrfinderplus.txt'
amrfinder_output        = 'amrfinderplus_mqc.txt'
blobtools_input         = 'blobtools_summary.txt'
blobtools_output        = 'blobtools_mqc.tsv'
circulocov_input        = 'circulocov_summary.tsv'
circulocov_output       = 'circulocov_mqc.tsv'
core_genome_input       = 'core_genome_evaluation.csv'
core_genome_output      = 'core_genome_evaluation_mqc.csv' 
drprg_input             = 'drprg_summary.tsv'
drprg_output            = 'drprg_mqc.tsv'
elgato_input            = 'elgato_summary.tsv'
elgato_output           = 'elgato_mqc.tsv'
emmtyper_input          = 'emmtyper_summary.tsv'
emmtyper_output         = 'emmtyper_mqc.tsv'
fastani_input           = 'fastani_summary.csv'
fastani_output          = 'fastani_mqc.csv'
heatcluster_input       = 'heatcluster.png'
heatcluster_output      = 'heatcluster_mqc.png'
kaptive_input           = 'kaptive_summary.txt'
kaptive_output          = 'kaptive_mqc.tsv'
kleborate_input         = 'kleborate_results.tsv'
kleborate_output        = 'kleborate_mqc.tsv'
mash_input              = 'mash_summary.csv'
mash_output             = 'mash_mqc.csv'
mlst_input              = 'mlst_summary.tsv'
mlst_output             = 'mlst_mqc.tsv'
mykrobe_input           = 'mykrobe_summary.csv'
mykrobe_output          = 'mykrobe_summary_mqc.csv'
pbp_input               = 'pbptyper_summary.tsv'
pbp_output              = 'pbptyper_mqc.tsv'
phytreeviz_iqtr_input   = 'iqtree_tree.png'
phytreeviz_iqtr_output  = 'phytreeviz_iqtree2_mqc.png' 
phytreeviz_mshtr_input  = 'mashtree_tree.png'
phytreeviz_mshtr_output = 'phytreeviz_mashtree_mqc.png'
plasmidfinder_input     = 'plasmidfinder_result.tsv'
plasmidfinder_output    = 'plasmidfinder_mqc.tsv' 
seqsero2_input          = 'seqsero2_results.txt'
seqsero2_output         = 'seqsero2_mqc.txt'
serotypefinder_input    = 'serotypefinder_results.txt'
serotypefinder_output   = 'serotypefinder_mqc.txt'
shigatyper_input        = 'shigatyper_summary.txt'
shigatyper_output       = 'shigatyper_mqc.txt'
shigatyper_hit_input    = 'shigatyper_hits.txt' 
shigatyper_hit_output   = 'shigatyper_hits_mqc.tsv'
snpdists_input          = 'snp_matrix.txt'
snpdists_output         = 'snpdists_matrix_mqc.txt' 

##########################################
# getting ready for multiqc              #
##########################################

if exists(blobtools_input) :
    blobtools_df = pd.read_table(blobtools_input)
    blobtools_df = blobtools_df[~blobtools_df['name'].isin(['all', 'no-hit', 'undef'])]

    samples = blobtools_df['sample'].drop_duplicates().tolist()
    organisms = sorted(blobtools_df['name'].drop_duplicates().tolist())
    blobtools_result_df = pd.DataFrame(columns=["sample"] + organisms)

    for sample in samples:
        result = [sample]
        for organism in organisms:
            readper = blobtools_df.loc[(blobtools_df['sample'] == sample) & (blobtools_df['name'] == organism), 'bam0_read_map_p']
            orgper = readper.iloc[0] if not readper.empty else 0
            result = result + [orgper]
        
        blobtools_result_df.loc[len(blobtools_result_df.index)] = result

    blobtools_result_df.to_csv(blobtools_output, index=False, sep="\t")

if exists(drprg_input):
    drprg_df = pd.read_table(drprg_input)
    drprg_df.to_csv(drprg_output, index=False, sep="\t")

if exists(circulocov_input):
    circulocov_df = pd.read_table(circulocov_input)
    circulocov_df = circulocov_df[circulocov_df['contigs'] == 'all']
    circulocov_df = circulocov_df.drop(['circ', 'contigs'], axis=1)
    circulocov_df.to_csv(circulocov_output, index=False, sep="\t")

if exists(elgato_input):
    elgato_df = pd.read_table(elgato_input)
    elgato_df.to_csv(elgato_output, index=False, sep="\t")

if exists(emmtyper_input):
    emmtyper_df = pd.read_table(emmtyper_input)
    emmtyper_df.to_csv(emmtyper_output, index=False, sep="\t")

if exists(mash_input) :
    mash_df = pd.read_csv(mash_input)
    samples = mash_df['sample'].drop_duplicates().tolist()
    organisms = sorted(mash_df['organism'].drop_duplicates().tolist())
    mash_result_df = pd.DataFrame(columns=["sample"] + organisms)

    for sample in samples:
        df_len = len(mash_result_df)
        mash_result_df.loc[df_len] = pd.Series()
        mash_result_df.at[df_len, "sample"] = sample
        sample_df = mash_df[mash_df['sample'] == sample].copy()

        for organism in organisms:
            mash_result_df.at[df_len, organism] = (sample_df['organism'] == organism).sum()
        #mash_result_df = mash_result_df.reset_index()

        # only worked in python 3.12 :(
        #sample_df = mash_df[mash_df['sample' ] == sample].copy()
        #counts_df = sample_df['organism'].value_counts().reset_index()
        #counts_df = counts_df.rename(columns={"index": "organism", 0: "count"}) 
        #counts_df = counts_df.set_index('organism')
        #counts_df.columns = sample
        #counts_df = counts_df.transpose()
        #mash_result_df = pd.concat([mash_result_df, counts_df], axis=0, join='outer')

    mash_result_df = mash_result_df.fillna(0)
    mash_result_df.to_csv(mash_output, index=False)

if exists(fastani_input):
    fastani_df = pd.read_csv(fastani_input)
    fastani_df.to_csv(fastani_output)

if exists(shigatyper_hit_input):
    shigatyper_df = pd.read_table(shigatyper_hit_input)
    shigatyper_df = shigatyper_df.drop(shigatyper_df.columns[1], axis=1)
    shigatyper_df = shigatyper_df.iloc[:, :2]
    shigatyper_df.to_csv(shigatyper_hit_output, sep="\t")
    #if [ -f 'shigatyper_results.txt' ]     ; then awk '{print \$1 "_" \$2 "\t" \$3}' shigatyper_results.txt > shigatyper_mqc.tsv ; fi

if exists(shigatyper_input):
    shigatyper_df = pd.read_table(shigatyper_input)
    shigatyper_df.to_csv(shigatyper_output, index=False, sep="\t")

if exists(amrfinder_input):
    amrfinder_df = pd.read_table(amrfinder_input)
    amrfinder_df = amrfinder_df.replace(' ', '_', regex=True)
    amrfinder_df.to_csv(amrfinder_output, sep="\t")

if exists(kaptive_input):
    kaptive_df = pd.read_table(kaptive_input)
    kaptive_df = kaptive_df.replace(' ', '_', regex=True)
    kaptive_df.to_csv(kaptive_output, index=True, sep="\t")


if exists(kleborate_input):
    kleborate_df = pd.read_table(kleborate_input)
    kleborate_df = kleborate_df.iloc[:, [1] + list(range(2, 12))]
    kleborate_df.to_csv(kleborate_output, index=False, sep="\t")
    #if [ -f 'kleborate_results.tsv' ]      ; then cut -f 1,3-12 kleborate_results.tsv > kleborate_mqc.tsv       ; fi

if exists(mlst_input):
    mlst_df = pd.read_table(mlst_input)
    mlst_df = mlst_df.replace(' ', '_', regex=True)
    mlst_df = mlst_df.loc[:, ['sample', 'matching PubMLST scheme', 'ST']] 
    mlst_df.to_csv(mlst_output, index=False, sep="\t")

if exists (pbp_input):
    pbp_df = pd.read_table(pbp_input)
    pbp_df.to_csv(pbp_output, index=False, sep="\t")

if exists(plasmidfinder_input):
    plasmidfinder_df = pd.read_table(plasmidfinder_input)
    plasmidfinder_df = plasmidfinder_df.iloc[:, :5]
    plasmidfinder_df.to_csv(plasmidfinder_output, sep="\t")
    #if [ -f 'plasmidfinder_result.tsv' ]   ; then awk '{ print \$1 "_" NR "\t" \$2 "\t" \$3 "\t" \$4 "\t" \$5 }' plasmidfinder_result.tsv > plasmidfinder_mqc.tsv   ; fi

if exists(seqsero2_input):
    seqsero2_df = pd.read_table(seqsero2_input)
    seqsero2_df = seqsero2_df.iloc[:, [0] + list(range(3, 10))]
    seqsero2_df.to_csv(seqsero2_output, index=False, sep="\t")
    #if [ -f 'seqsero2_results.txt' ]       ; then cut -f 1,4-10 seqsero2_results.txt > seqsero2_mqc.txt ; fi

if exists(serotypefinder_input):
    serotypefinder_df = pd.read_table(serotypefinder_input)
    serotypefinder_df = serotypefinder_df.iloc[:, :6]
    serotypefinder_df = serotypefinder_df.replace(to_replace=' ', value='', regex=True)
    serotypefinder_df.columns = serotypefinder_df.columns.str.replace(' ', '')
    serotypefinder_df.to_csv(serotypefinder_output, index=True, sep="\t")
    #if [ -f 'serotypefinder_results.txt' ] ; then cut -f 1-6 serotypefinder_results.txt > serotypefinder_mqc.txt  ; fi
   
if exists(core_genome_input):
    core_genome_df = pd.read_csv(core_genome_input)
    core_genome_df = core_genome_df.loc[:, ['sample', 'core', 'soft', 'shell', 'cloud']]
    core_genome_df.to_csv(core_genome_output, index=False)

if exists(heatcluster_input):
    shutil.copyfile(heatcluster_input, heatcluster_output)

if exists(snpdists_input):
    shutil.copyfile(snpdists_input, snpdists_output)

if exists(phytreeviz_iqtr_input):
    shutil.copyfile(phytreeviz_iqtr_input, phytreeviz_iqtr_output)

if exists(phytreeviz_mshtr_input):
    shutil.copyfile(phytreeviz_mshtr_input, phytreeviz_mshtr_output)

if exists(mykrobe_input):
    shutil.copyfile(mykrobe_input, mykrobe_output)
