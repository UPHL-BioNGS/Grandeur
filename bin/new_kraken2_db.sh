#/bin/bash

#################################################
# written by Erin Young                         #
# for downloading kraken2 database for grandeur #
#################################################

############################################################

USAGE="
Downloads Standard 8 Kraken2 database (https://benlangmead.github.io/aws-indexes/k2)
with wget.

Usage: new_kraken2_db.sh
"

############################################################

echo "$USAGE"

mkdir kraken2_db
cd kraken2_db

echo "$(date): Downloading Kraken2 database"
wget -c "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20221209.tar.gz"

echo "$(date): Extracting Kraken2 database"

tar -zxvf k2_standard_08gb_20221209.tgz
echo "$(date): Finished"

############################################################