#/bin/bash

# just a bunch of tests with local directories
# /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur/bin/.tests.sh

# default with input
nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur \
  -profile singularity \
  --sample_sheet /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur/bin/sample_sheet.csv \
  --outdir grandeur_sample_sheet \
  -resume  \
  -with-tower

# default with reads
nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur \
  -profile singularity \
  --reads  /home/eriny/sandbox/test_files/grandeur/reads \
  --outdir grandeur_fastq_channel \
  -resume  \
  -with-tower

# default with fastas
nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur \
  -profile singularity \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir grandeur_fasta_channel \
  -resume  \
  -with-tower

# default with reads and fastas
nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur \
  -profile singularity \
  --reads  /home/eriny/sandbox/test_files/grandeur/reads \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir grandeur_fastq_fasta_channel \
  -resume  \
  -with-tower

# multiple sequence alignment
nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur \
  -profile singularity,msa \
  --gff    /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --reads  /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir grandeur_msa_fastani \
  --min_core_genes 50 \
  -resume  \
  -with-tower

for profile in "singularity" "uphl"
do
  for ver in "test0" "test1" "test2" "test3" "test4" "test5" "test6"
  do
    nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur \
      -profile $profile,$ver \
      --outdir grandeur_${ver}_$profile \
      -resume  \
      -with-tower
  done
done

# with nothing
nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Grandeur \
  -profile singularity \
  --gff    wontexist \
  --fastas shouldntexit \
  --reads  doesntexist \
  --outdir grandeur_empty \
  -resume  \
  -with-tower
