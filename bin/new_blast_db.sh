#/bin/bash

#################################################
# written by Erin Young                         #
# for downloading kraken2 database for grandeur #
#################################################

############################################################

USAGE="
Downloads RefSeq database for blast
with wget.

Usage: 
new_blast_db.sh
new_blast_db.sh destination_dir
"

############################################################

# Dependency check
if [ -z "$(which wget)" ] ; then echo "$(date) : FATAL : wget not found!" ; exit 1 ; fi
if [ -z "$(which tar)" ]  ; then echo "$(date) : FATAL : tar not found!"  ; exit 1 ; fi

############################################################

# Setting final directory
if [ -z "$1" ] ; then out=blast_db ; else out=$1 ; fi

# creating the directory
echo "$(date): Starting blast download in $out"
mkdir -p $out
cd $out

# getting taxbd
wget --continue --show-progress "https://ftp.ncbi.nlm.nih.gov/blast/db/v5/taxdb.tar.gz"

echo "$(date): Downloading tar files from NCBI"

# yes, there are too many i's - it's so I don't have to update later
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
do
    wget --continue --show-progress "https://ftp.ncbi.nlm.nih.gov/blast/db/v5/ref_prok_rep_genomes.$i.tar.gz" || break
    tar -xvf ref_prok_rep_genomes.$i.tar.gz
    rm ref_prok_rep_genomes.$i.tar.gz
done

echo "$(date): File preparation is complete. Blast database is in $out."

# works on mac but not LW
# wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_prok_rep_genomes.??.tar.gz"

############################################################
