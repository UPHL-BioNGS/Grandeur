#!/usr/bin/env python3

'''
Author: Erin Young

Description:

This script is to get some genome accession from NCBI datasets

EXAMPLE:
python3 evaluate.py
'''

import pandas as pd
from pathlib import Path

genepresence_df = pd.read_table("gene_presence_absence.Rtab")
num_samples = len(genepresence_df.columns) - 1

genepresence_df['sum'] = genepresence_df.drop('Gene', axis=1).sum(axis=1)

core_df = genepresence_df[ genepresence_df['sum'] >= num_samples * .99 ]
soft_df = genepresence_df[(genepresence_df['sum'] >= num_samples * .95 ) & (genepresence_df['sum'] < num_samples * .99 )]
shel_df = genepresence_df[(genepresence_df['sum'] >= num_samples * .15 ) & (genepresence_df['sum'] < num_samples * .95 )]
clud_df = genepresence_df[ genepresence_df['sum'] <  num_samples * .15 ]

samples = genepresence_df.drop('Gene', axis=1).drop('sum', axis=1).columns
percent_df = pd.DataFrame([])
percent_df['sample'] = samples

for sample in samples:
    bamindex = percent_df.index[percent_df['sample'] == sample]
    total = genepresence_df[sample].sum(axis=0)
    core  = core_df[sample].sum(axis=0)
    soft  = soft_df[sample].sum(axis=0)
    shell = shel_df[sample].sum(axis=0)
    cloud = clud_df[sample].sum(axis=0)
    
    percent_df.loc[bamindex, 'total'] = total
    percent_df.loc[bamindex, 'core']  = core
    percent_df.loc[bamindex, 'soft']  = soft
    percent_df.loc[bamindex, 'shell'] = shell
    percent_df.loc[bamindex, 'cloud'] = cloud


percent_df["per_core"]  = percent_df['core']  / percent_df['total']
percent_df["per_soft"]  = percent_df['soft']  / percent_df['total']
percent_df["per_shell"] = percent_df['shell'] / percent_df['total']
percent_df["per_clouc"] = percent_df['cloud'] / percent_df['total']
percent_df = percent_df.sort_values('per_core', ascending=False)

core_genome_file = Path("core_gene_alignment_filtered.aln")
if not core_genome_file.is_file():
    core_genome_file = Path("core_gene_alignment.aln")

sample    = ""
length    = ""
ambiguous = ""
with open(core_genome_file) as file:
    for line in file:
        if ">" in line:
            if sample:
                bamindex = percent_df.index[percent_df['sample'] == sample]
                percent_df.loc[bamindex, 'length']        = length
                percent_df.loc[bamindex, 'num_ambiguous'] = ambiguous               
            sample     = line.replace(">","").strip()
            ambiguous  = 0
            length     = 0
        else:
            line       = line.strip()
            length    += len(line)
            nonagct    = len(line) - line.count('a') - line.count('A') - line.count('g') - line.count('G') - line.count('c') - line.count('C') - line.count('t') - line.count('T')
            ambiguous += nonagct

bamindex = percent_df.index[percent_df['sample'] == sample]
percent_df.loc[bamindex, 'length']        = length
percent_df.loc[bamindex, 'num_ambiguous'] = ambiguous

percent_df["per_ambiguous"] = percent_df['num_ambiguous'] / percent_df['length']

percent_df.to_csv('core_genome_evaluation.csv', index=False)
