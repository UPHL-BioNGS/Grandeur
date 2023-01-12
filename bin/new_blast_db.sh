#/bin/bash

#################################################
# written by Erin Young                         #
# for downloading kraken2 database for grandeur #
#################################################

############################################################

USAGE="
Downloads RefSeq database for blast
with wget.

Usage: new_blast_db.sh
"

############################################################

echo "$(date): Starting blast download"
mkdir -p blast_db
cd blast_db

# wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"

wget --continue --show-progress "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/taxdb.tar.gz"

# works on mac but not LW
# wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_prok_rep_genomes.??.tar.gz"

echo "$(date): Downloading tar files from NCBI"

# yes, there are too many i's
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
do
    wget --continue --show-progress "https://ftp.ncbi.nlm.nih.gov/blast/db/v5/ref_prok_rep_genomes.$i.tar.gz"
    tar -xvf ref_prok_rep_genomes.$i.tar.gz
    rm ref_prok_rep_genomes.$i.tar.gz
done

# too much for the LW network
# i="00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20"
#
# echo $i | tr " " "\n" | parallel "wget --continue --show-progress \"https://ftp.ncbi.nlm.nih.gov/blast/db/v5/ref_prok_rep_genomes.{}.tar.gz\""
#
# echo "$(date): Extracting files and removing tar.gz files"
# ls *tar.gz | parallel "tar -xvf {} && rm {} "

echo "$(date): File preparation is complete"