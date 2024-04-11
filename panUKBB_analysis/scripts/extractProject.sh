#!/bin/bash
#SBATCH -p shared
#SBATCH -J GLEANER_full_test
#SBATCH --time=60:00
#SBATCH --mem=50G
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization
ml anaconda
conda activate std

python src/quickGWASIter.py  --type ldsc_custom --output gwas_extracts/panUKBB_complete/full_hapmap3_snps --gwas_list gwas_extracts/panUKBB_complete/missingness_report.tsv --snp_list /data/abattle4/aomdahl1/reference_data/ldsc_tutorial_hm3.no_hla.snplist --extension .sumstats.gz --gwas_dir  gwas_extracts/panUKBB_complete/
#How muchg missingness here
#Project the sumstats, with correction (arguments updated)
conda deactivate
conda activate renv


#Project the sumstats, with correction (arguments updated)
python src/projectSumStats.R --factors /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_6.1K/latent.factors.txt --sumstats /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/ull_hapmap3_snps.z.tsv --id_type RSID --output "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_6.1K/projectedMetaAdj_hapmap3_loadings.txt"  --proj_method meta --trait_names "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete_6.1K/trait_out_order.txt" --decorrelate /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_6.1K/_scaledShrunkBlockCovarMatrix.txt

#Proect sumstats, without correction:
python src/projectSumStats.R --factors /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_6.1K/latent.factors.txt --sumstats /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/hm3_outputull_hapmap3_snps.z.tsv --id_type RSID --output "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_6.1K/projectedMetaAdj_hapmap3_loadings.txt"  --proj_method meta --trait_names "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete_6.1K/trait_out_order.txt"

#Move some files around to make sure they are in teh right spot for snakemake
mkdir -p gwas_extracts/panUKBB_complete_6.1K/
cp gwas_extracts/panUKBB_complete/full_hapmap3_snps* gwas_extracts/panUKBB_complete_6.1K/
