#/bin/bash

# the default workflow on just fastas
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity,just_fasta \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir default_fastas

# the default workflow on just paired-end reads
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity,just_fastq \
  --reads /home/eriny/sandbox/test_files/grandeur/reads \
  --outdir default_reads \
  -resume

# the default workflow on paired-end reads plus fastas
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --reads /home/eriny/sandbox/test_files/grandeur/reads \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir default_reads_fastas \
  -resume

# the default workflow on paired-end reads plus fastas with fastq_to_consensus profile
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity,fastq_to_consensus \
  --reads /home/eriny/sandbox/test_files/grandeur/reads \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir default_reads_profile_fastas \
  -resume

# the default workflow on paired-end reads plus fastas with uphl's config
nextflow run /home/eriny/sandbox/Grandeur \
  -profile uphl \
  --reads /home/eriny/sandbox/test_files/grandeur/reads \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir default_uphl \
  -resume

# doing a multiple sequence alignment with gff, fastqs, and fastas
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity,msa \
  --gff /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --reads /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir everything \
  --prokka true \
  --roary true \
  -resume

# doing a multiple sequence alignment with gff, fastqs, and fastas with fastq_to_msa profile
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity,fastq_to_msa \
  --gff /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --reads /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir everything_profile \
  -resume

# doing a multiple sequence alignment with gff, fastqs, and fastas with fastq_to_msa and extras_off profiles
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity,fastq_to_msa,extras_off \
  --gff /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --reads /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir everything_profile_2 \
  -resume

# just doing a multiple sequence alignment with gff and fastas
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --gff /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir msa \
  --prokka true \
  --roary true \
  -resume

# just doing a multiple sequence alignment with gff and fastas with fasta_to_msa profile
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity,fasta_to_msa \
  --gff /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir msa_profile \
  -resume

# with nothing
nextflow run /home/eriny/sandbox/Grandeur \
  -profile singularity \
  --gff wontexist \
  --fastas shouldntexit \
  --reads doesntexist \
  --outdir empty \
  -resume
