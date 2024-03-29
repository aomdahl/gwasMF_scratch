#!/bin/bash
#SBATCH -p express
#SBATCH -J second_half_extraction
#SBATCH --time=90:00
source /data/apps/go.sh ### for safety reasons
set -e
cd /data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files
bash extract_EUR_scripts/phecode-272.11-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-272.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-300-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-317-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-351-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-361.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-361-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-365-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-366-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-401.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-401-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-411.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-411.2-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-411.4-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-411-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-427.2-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-427-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-454.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-454-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-455-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-471-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-495-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-496-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-530.11-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-530.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-530-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-531-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-550.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-550.4-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-550-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-562.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-562-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-574.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-574-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-594-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-726-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-728.71-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-728.7-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-735.3-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-740.1-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-740-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-788-both_sexes.tsv.extract.sh
bash extract_EUR_scripts/phecode-835-both_sexes.tsv.extract.sh
