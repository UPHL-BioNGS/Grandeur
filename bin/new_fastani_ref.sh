#/bin/bash

############################################
# written by Erin Young                    #
# for creating reference file for grandeur #
############################################

# Setting final directory
if [ -z "$1" ] ; then out=genomes ; else out=$1 ; fi

if [ -d "$out" ] ; then echo "$out directory already exists" ; exit 1 ; fi

if [ -z "$(which datasets)" ] ; then echo "$(date) : FATAL : datasets not found!" ; exit 1 ; fi

mkdir $out
cd $out

accessions="GCF_008632635.1 GCF_000191145.1 GCF_001558935.2 GCF_003812345.1 GCF_003812345.1 GCF_000018045.1 GCF_900638065.1 GCF_002023665.2 GCF_007035805.1 GCF_015137655.1 GCF_023702375.1 GCF_023702375.1 GCF_019048625.1 GCF_000534275.1 GCA_002741475.1 GCF_007632255.1 GCF_015139575.1 GCF_003812925.1 GCF_000240185.1 GCF_016415705.1 GCF_009648975.1 GCF_902387845.1 GCF_003019925.1 GCF_000069965.1 GCF_003204135.1 GCF_023547145.1 GCF_000006765.1 GCF_901421005.1 GCF_003516165.1 GCF_900475405.1"

datasets download genome accession $accessions --filename ncbi_dataset.zip

# cut -f 1 /home/eriny/sandbox/Grandeur/configs/genomes.txt > id_list.txt
datasets download genome accession --inputfile id_list.txt --filename ncbi_dataset.zip

unzip -o ncbi_dataset.zip
    
fastas=$(ls ncbi_dataset/data/*/*.fna)

for fasta in ${fastas[@]}
do  
    accession=$(echo $fasta | cut -f 3 -d / )
    organism=$(head -n 1 $fasta | awk '{print $2 "_" $3 }' )
    cat $fasta | sed 's/ /_/g' | sed 's/,//g' > $out/${organism}_${accession}.fna
done

rm ncbi_dataset.zip
rm -rf ncbi_dataset
cd ../

tar -czvf fastani_refs.tar.gz $out/

# mv fastani_refs.tar.gz ~/sandbox/Grandeur/configs/fastani_refs.tar.gz