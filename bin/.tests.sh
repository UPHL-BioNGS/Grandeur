#/bin/bash

# just a bunch of tests with local directories
# /home/eriny/sandbox/Grandeur/bin/.tests.sh

# default with reads
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --reads  /home/eriny/sandbox/test_files/grandeur/reads \
  --outdir fastq_channel \
  -resume  \
  -with-tower

# default with fastas
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir fasta_channel \
  -resume  \
  -with-tower

# default with reads and fastas
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --reads  /home/eriny/sandbox/test_files/grandeur/reads \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir fastq_fasta_channel \
  -resume  \
  -with-tower

# multiple sequence alignment
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity,msa \
  --gff    /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --reads  /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir msa \
  -resume  \
  -with-tower

for profile in "singularity" "uphl"
do
  # test
  nextflow run /home/eriny/sandbox/Grandeur \
    -profile $profile,test \
    --outdir test_$profile \
    -resume  \
    -with-tower

  # test1
  nextflow run /home/eriny/sandbox/Grandeur \
    -profile $profile,test1 \
    --outdir test1_$profile \
    -resume  \
    -with-tower

  # test2
  nextflow run /home/eriny/sandbox/Grandeur \
    -profile $profile,test2 \
    --outdir test2_$profile \
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
