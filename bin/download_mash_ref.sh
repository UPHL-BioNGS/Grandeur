#/bin/bash

################################################
# written by Erin Young                        #
# for downloading Erin's msh file for grandeur #
# refseq 217                                   #
################################################

ver=215
URL=https://zenodo.org/record/7887021/files/rep-genomes.msh

############################################################

USAGE="
Downloads a premade mash reference file for refseq v$ver
from $URL

download_mash_ref.sh mash_v$ver
"

############################################################

echo "$USAGE"

# Setting final directory
if [ -z "$1" ] ; then out=mash_v$ver ; else out=$1 ; fi

echo "$(date): Downloading file"

mkdir $out
cd $out

wget --continue --show-progress $URL

echo "$(date): Finished!"

############################################################
