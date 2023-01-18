#/bin/bash

#################################################
# written by Erin Young                         #
# for downloading kraken2 database for grandeur #
#################################################

VER="k2_standard_08gb_20221209"

############################################################

USAGE="
Downloads $VER pre-built Kraken2 database, a standard 8G Kraken2 database, 
with wget and extract with tar.

(More databases are listed at https://benlangmead.github.io/aws-indexes/k2)

Usage: new_kraken2_db.sh kraken2_db
"

############################################################

echo "$USAGE"

############################################################

# Setting final directory
if [ -z "$1" ] ; then out=kraken2_db ; else out=$1 ; fi

# creating the directory
echo "$(date): Kraken2 database will be in $out"
mkdir -p $out
cd $out

echo "$(date): Downloading Kraken2 database"
wget -c "https://genome-idx.s3.amazonaws.com/kraken/$VER.tar.gz"

echo "$(date): Extracting Kraken2 database"

tar -zxvf $VER.tar.gz || echo "$(date) : WARNING : could not extract $VER.tar.gz "
echo "$(date): Finished"

############################################################