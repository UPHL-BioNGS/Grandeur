#/bin/bash

################################################
# written by Erin Young                        #
# for downloading Erin's msh file for grandeur #
# refseq 215                                   #
################################################

ver=215

############################################################

USAGE="
Downloads a premade mash reference file for refseq v$ver
"

############################################################

echo "$USAGE"

echo "$(date): Downloading file"

mkdir mash_v$ver
cd mash_v$ver

wget --continue --show-progress $URL

echo "$(date): Finished!"

############################################################
