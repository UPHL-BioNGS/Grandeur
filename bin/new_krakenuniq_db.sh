#/bin/bash

####################################################
# written by Erin Young                            #
# for downloading krakenuniq database for grandeur #
####################################################

############################################################

USAGE="
Downloads ReSeq for krakenuniq database with krakenuniq-download.

Usage: new_krakenuniq_db.sh
"

############################################################

echo "$USAGE"


if [ -z "$(which krakenuniq-download)" ] ; then echo "$(date): krakenuniq-download not found!" ; exit 1 ; fi

mkdir -p krakenuniq_db

echo "$(date): downloading krakenuniq database"

krakenuniq-download --db krakenuniq_db --threads 10 --dust refseq/bacteria

echo "$(date): building krakenuniq database"

# The building step may take up to a couple of days on large sequence sets such as nt.

krakenuniq-build --db krakenuniq_db --kmer-len 31 --threads 10 --taxids-for-genomes --taxids-for-sequences

echo "$(date): finished!"