#!/bin/bash
#$ -cwd # Run job from current directory
#$ -N grandeur_refactor # Job name
#$ -l h_vmem=128G # Memory requirement
#$ -l h_rt=04:00:00 # Max runtime (4 hours)
#$ -pe smp 8 # Request 8 CPUs
#$ -M xre6@cdc.gov
#$ -m abe
#$ -q short.q

source /etc/profile

module purge
module load nextflow
module load singularity

nextflow run /scicomp/home-pure/xre6/Projects/grandeur/main.nf \
-profile my_test2,singularity \
--outdir /scicomp/home-pure/xre6/Projects/grandeur/grandeur_testing/read_identification_refactor \
--blast_db /scicomp/reference-pure/ncbi-blast-databases/ \
--blast_db_type core_nt
