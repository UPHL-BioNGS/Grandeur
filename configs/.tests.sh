#/bin/bash

# one channel at a time
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --reads  /home/eriny/sandbox/test_files/grandeur/reads \
  --outdir fastq_channel \
  -resume  \
  -with-tower

nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir fasta_channel \
  -resume  \
  -with-tower

nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --gff    /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir gff_channel \
  -resume  \
  -with-tower

# testing profiles
for profile in "fastq_to_consensus" "just_fastq" "just_fasta" "extras_off" "uphl"
do
  nextflow run /home/eriny/sandbox/Grandeur \
    -profile singularity,$profile \
    --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
    --reads  /home/eriny/sandbox/test_files/grandeur/reads \
    --outdir profile_$profile \
    -resume  \
    -with-tower
done

# testing msa-related profiles
for profile in "msa" "extras_off,msa" "fastq_to_msa" "fasta_to_msa"
do
  nextflow run /home/eriny/sandbox/Grandeur \
    -profile singularity,$profile \
    --gff    /home/eriny/sandbox/test_files/grandeur/msa \
    --fastas /home/eriny/sandbox/test_files/grandeur/msa \
    --reads  /home/eriny/sandbox/test_files/grandeur/msa \
    --outdir msa_$profile \
    -resume  \
    -with-tower
done

# with nothing
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --gff    wontexist \
  --fastas shouldntexit \
  --reads  doesntexist \
  --outdir empty \
  -resume  \
  -with-tower
