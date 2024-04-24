#!/bin/bash
#SBATCH -p shared
#SBATCH -J simsnoLEANER_full_test
#SBATCH --time=180:00
#SBATCH --mem=15G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style
ml anaconda
conda activate renv
#Rscript gwas_centric_non_null_zscores_fullDat.R 3 1
Rscript gwas_centric_non_null_zscores_100.R   4  1
conda deactivate
