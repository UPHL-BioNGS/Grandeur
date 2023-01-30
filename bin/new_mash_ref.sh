#/bin/bash

######################################
# written by Erin Young              #
# for creating msh file for grandeur #
######################################

############################################################

USAGE="
Uses NCBI's datasets and dataformat to download 
representative genomes for bacteria, combine those sequences 
into one fasta, and then create a mash sketch.

Usage: new_mash_ref.sh
"

############################################################

echo "$USAGE"

if [ -z "$(which datasets)" ]   ; then echo "$(date): FATAL : datasets not found"   ; exit 1 ; fi
if [ -z "$(which dataformat)" ] ; then echo "$(date): FATAL : dataformat not found" ; exit 1 ; fi
if [ -z "$(which mash)" ]       ; then echo "$(date): FATAL : mash not found"       ; exit 1 ; fi

# Setting final directory
if [ -z "$1" ] ; then out=new_mash ; else out=$1 ; fi

echo "$(date): Getting ids for representative genomes"

mkdir $out
cd $out

datasets summary genome taxon bacteria --reference --as-json-lines | \
  dataformat tsv genome --fields accession,assminfo-refseq-category,organism-name --elide-header | \
  grep representative | \
  tee representative_genomes.txt | \
  cut -f 1 > genome_ids.txt

echo "$(date): Downloading genomes for ids"
datasets download genome accession --inputfile genome_ids.txt --filename rep-genomes.zip

echo "$(date): Decompressing zip file"
unzip rep-genomes.zip

echo "$(date): Creating file for mash"
cat  ncbi_dataset/data/*/*.fna  | sed 's/ /_/g' | sed 's/,//g' > rep-genomes.fasta

echo "$(date): Skeching rep-genomes.fasta"
mash sketch -i -p 20 rep-genomes.fasta -o rep-genomes

############################################################

echo "$(date): File preparation is complete"
ls -alh rep-genomes.fasta
ls -alh rep-genomes.msh

echo "$(date): Remaining tasks:
- upload rep-genomes.msh to Zenodo
- update the github readme
"

############################################################
