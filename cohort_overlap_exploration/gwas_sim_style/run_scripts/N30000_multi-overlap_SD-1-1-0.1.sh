#!/bin/bash
#SBATCH -p shared
#SBATCH -J process_sims 
#SBATCH --time=10:00:00
#SBATCH --mem=30G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style
ml anaconda
conda activate renv
set -e
mkdir -p /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/N30000_multi-overlap_SD-1-1-0.1_h2-0-0/
echo "Generating simulations"
Rscript simulation_scripts/gwas_centric_vary_overlap_degree_template.R 6 10 N30000_multi-overlap_SD-1-1-0.1_h2-0-0/  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/ 0 0 1,1,0.1

echo "Simulation data generated. Evaluating now..."
Rscript simulation_scripts/factorize_input_simulations.R --input_dir /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/N30000_multi-overlap_SD-1-1-0.1_h2-0-0/ --output /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/N30000_multi-overlap_SD-1-1-0.1_h2-0-0/factorization.summaries.RData --with_covar
conda deactivate
