#!/bin/bash
#SBATCH -p shared
#SBATCH -J process_sims 
#SBATCH --time=180:00
#SBATCH --mem=10G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style
ml anaconda
conda activate renv
#Rscript simulation_scripts/factorize_input_simulations.R --input_dir /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/large_studies_with_covar/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/large_studies_with_covar/factorization.summaries.RData  --with_covar
#Rscript simulation_scripts/factorize_input_simulations.R --input_dir /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/large_studies_with_covar/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/large_studies_with_covar/factorization.summaries.post750.RData  --with_covar --start 751
Rscript simulation_scripts/factorize_input_simulations.R --input_dir /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/large_studies_with_covar/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/large_studies_with_covar/factorization.summaries.post1750.RData  --with_covar --start 1751
conda deactivate
