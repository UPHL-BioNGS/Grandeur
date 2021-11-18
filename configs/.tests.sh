#/bin/bash

# the default workflow on just paired-end reads
nextflow /home/eriny/sandbox/Grandeur/grandeur.nf \
  -profile singularity \
  --reads /home/eriny/sandbox/test_files/grandeur/reads \
  --outdir default_reads \
  -with-tower --shigatyper false

# the default workflow on just fastas
nextflow /home/eriny/sandbox/Grandeur/grandeur.nf \
  -profile singularity \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir default_fastas \
  -with-tower --shigatyper false

# the default workflow on paired-end reads plus fastas
nextflow /home/eriny/sandbox/Grandeur/grandeur.nf \
  -profile singularity \
  --reads /home/eriny/sandbox/test_files/grandeur/reads \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  --outdir default_reads_fastas \
  -with-tower --shigatyper false

# the default workflow on paired-end reads plus fastas with uphl's config
nextflow /home/eriny/sandbox/Grandeur/grandeur.nf -profile singularity,uphl \
  --reads /home/eriny/sandbox/test_files/grandeur/reads \
  --fastas /home/eriny/sandbox/test_files/grandeur/fastas \
  -resume \
  --outdir default_uphl \
  -with-tower --shigatyper false

# doing a multiple sequence alignment with gff, fastqs, and fastas
nextflow /home/eriny/sandbox/Grandeur/grandeur.nf \
  -profile singularity \
  --gff /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --reads /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir everything \
  --roary true \
  -resume \
  -with-tower

# just doing a multiple sequence alignment with gff and fastas
nextflow /home/eriny/sandbox/Grandeur/grandeur.nf \
  -profile singularity,msa \
  --gff /home/eriny/sandbox/test_files/grandeur/msa \
  --fastas /home/eriny/sandbox/test_files/grandeur/msa \
  --outdir msa \
  -resume \
  -with-tower

# with nothing
nextflow /home/eriny/sandbox/Grandeur/grandeur.nf \
  -profile singularity \
  --gff wontexist \
  --fastas shouldntexit \
  --reads doesntexist \
  --outdir empty \
  -with-tower
