#!/bin/bash
#SBATCH -p shared
#SBATCH -J simsquick
#SBATCH --time=240:00
#SBATCH --mem=25G

source /data/apps/go.sh ###

cd /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style
ml anaconda
conda activate renv
Rscript simulation_scripts/gwas_centric_vary_overlap_degree_template.R 6 50 large_studies_with_covar/  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/ 0 0

conda deactivate
