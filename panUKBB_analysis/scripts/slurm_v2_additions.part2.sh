#!/bin/bash
#SBATCH -p express
#SBATCH -J first_half_extraction
#SBATCH --time=90:00
source /data/apps/go.sh ### for safety reasons
set -e
cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files
bash extract_EUR_scripts/categorical-20110-both_sexes-9.tsv.extract.sh
bash extract_EUR_scripts/categorical-20117-both_sexes-1.tsv.extract.sh
bash extract_EUR_scripts/categorical-2492-both_sexes-2492.tsv.extract.sh
bash extract_EUR_scripts/categorical-4598-both_sexes-4598.tsv.extract.sh
bash extract_EUR_scripts/categorical-6159-both_sexes-4.tsv.extract.sh
bash extract_EUR_scripts/continuous-20519-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/continuous-22677-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-30150-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-4548-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/continuous-4825-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/continuous-51-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-5262-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-AG-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-DBP-both_sexes-combined_medadj_irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-eGFRcys-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-FEV1FVC-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-IBil-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-LDLC-both_sexes-medadj_irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-MAP-both_sexes-combined_irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-MCP-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-NAP-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-PP-both_sexes-combined_medadj_irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-SBP-both_sexes-combined_irnt.tsv.extract.sh
bash extract_EUR_scripts/continuous-whr-both_sexes-irnt.tsv.extract.sh
bash extract_EUR_scripts/icd10-C34-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/icd10-H40-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/icd10-I21-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/icd10-K21-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-172-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-174-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-208-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-250.2-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-250-both_sexes.tsv.extract.sh
