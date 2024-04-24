#!/bin/bash
#SBATCH -p shared
#SBATCH -J simsquick
#SBATCH --time=240:00
#SBATCH --mem=15G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style
ml anaconda
conda activate renv
Rscript simulation_scripts/gwas_centric_lower2_herit_zscores_2024.R 6 1 50
conda deactivate
