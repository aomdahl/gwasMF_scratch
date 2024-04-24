#!/bin/bash
#SBATCH --job-name=gwas_centric_sim
#SBATCH -p defq
#SBATCH -n 1
#SBATCH --time=3:10:0
#set -euo  pipefail
ml anaconda
conda activate renv
ml gsl
Rscript gwas_centric_spearman.R
conda deactivate
