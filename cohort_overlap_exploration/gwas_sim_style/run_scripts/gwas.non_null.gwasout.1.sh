#!/bin/bash
#SBATCH --job-name=gwas_centric_si1m
#SBATCH -p defq
#SBATCH -n 1
#SBATCH --time=3:00:0
#set -euo  pipefail
ml anaconda
conda activate renv
ml gsl
Rscript gwas_centric_non_null_zscores_100.R 1
conda deactivate
