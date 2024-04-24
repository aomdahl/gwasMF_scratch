#!/bin/bash
#SBATCH -p shared
#SBATCH -J process_sims 
#SBATCH --time=10:00:00
#SBATCH --mem=15G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style
ml anaconda
conda activate renv
#Rscript simulation_scripts/factorize_input_simulations.R --input_dir /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/large_studies_with_covar/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/large_studies_with_covarfactorization.summaries.RData  --with_covar
Rscript simulation_scripts/factorize_input_simulations.R --input_dir /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/large_studies_with_covar/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/large_studies_with_covarfactorization.summaries.RData  --with_covar --start 51
conda deactivate
