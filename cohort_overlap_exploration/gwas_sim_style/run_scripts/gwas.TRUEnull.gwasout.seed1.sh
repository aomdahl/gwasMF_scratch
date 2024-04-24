#!/bin/bash
#SBATCH --job-name=gwas_centric_sim_truenull
#SBATCH -p defq
#SBATCH -n 1
#SBATCH --time=3:00:0
#set -euo  pipefail
ml anaconda
conda activate renv
ml gsl
Rscript gwas_centric_TRUEnull_zscores_100.R 1
conda deactivate
