#!/bin/bash
#SBATCH -p shared
#SBATCH -J simsquick
#SBATCH --time=240:00
#SBATCH --mem=25G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style
ml anaconda
conda activate renv
#Rscript simulation_scripts/gwas_centric_lower2_herit_zscores_2024.R 6 1 50
Rscript simulation_scripts/gwas_centric_vary_overlap_degree_template.R 5 50 large_studies_with_covar  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/ 0.05 0.05

conda deactivate
