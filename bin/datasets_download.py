#!/usr/bin/env python3

'''
Author: Erin Young

Description:

This script is to get some genome accession from NCBI datasets

EXAMPLE:
python3 datasets_download.py taxon hits
'''

import subprocess
import sys

taxon = sys.argv[1]
genus, species = taxon.replace('[', '').replace(']', '').split('_')
print("Looking for accessions for " + genus + " " + species )
outfile = open('datasets/' + genus + "_" + species + '_genomes.csv', "w")

try:
  hits = sys.argv[2]
except:
  hits = '5'

# putting in the header
outfile.write('accession,assminfo-refseq-category,assminfo-level,organism-name,assmstats-total-ungapped-len\n')

# Getting representative genomes
rep = subprocess.Popen(['datasets', 'summary', 'genome', 'taxon', '"' + genus + ' ' + species + '"', '--reference', '--limit', hits, '--as-json-lines'], stdout = subprocess.PIPE)
dft = subprocess.check_output(['dataformat', 'tsv', 'genome', '--fields' , 'accession,assminfo-refseq-category,assminfo-level,organism-name,assmstats-total-ungapped-len'], stdin = rep.stdout, universal_newlines= True, text='str')
rep.wait()
for line in dft.split('\n'):
  if 'Ungapped Length' not in line and line:
    if int(line.split('\t')[4]) < 15000000:
      outfile.write(line.replace('\t',',') + '\n')

# Getting additional genomes
oth = subprocess.Popen(['datasets', 'summary', 'genome', 'taxon', '"' + genus + ' ' + species + '"', '--limit', hits, '--as-json-lines'], stdout = subprocess.PIPE)
df2 = subprocess.check_output(['dataformat', 'tsv', 'genome', '--fields' , 'accession,assminfo-refseq-category,assminfo-level,organism-name,assmstats-total-ungapped-len'], stdin = oth.stdout, universal_newlines= True, text='str')
oth.wait()
for line in df2.split('\n'):
  if 'Ungapped Length' not in line and line:
    if int(line.split('\t')[4]) < 15000000:
      outfile.write(line.replace('\t',',') + '\n')

outfile.close()
