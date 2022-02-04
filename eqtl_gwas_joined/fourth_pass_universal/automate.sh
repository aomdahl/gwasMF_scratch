#!/bin/bash
#Arguments are:
    # 1 SNP list
    # ....and that's it
SNPS=$1 #SNP list
#SNPS=./gwas_snp_lists/gwas_centric_1e-5.ids.txt
I=`basename $SNPS`
H="${I%.*}" #remove .txt
DIR=`echo $SNPS | awk -F "/" '{NF--; print}' | tr " " "/"`
set -e
GD=/work-zfs/abattle4/ashton/snp_networks/scratch/eqtl_gwas_joined/src/
echo $GD
echo $I
echo $H
#SRC=$2 #specify its its gwas, rand, or sq
#Prep
ml gcc/5.5.0
ml R
ml python/3.7-anaconda
#Clean up the list
echo "Cleaning up the list"
mkdir -p init_snp_lists
echo "Rscript /work-zfs/abattle4/ashton/snp_utils/ambig_multi_indel_cleanup.R $SNPS init_snp_lists/${H}.cleaned.txt"
#Rscript /work-zfs/abattle4/ashton/snp_utils/ambig_multi_indel_cleanup.R $SNPS init_snp_lists/${H}.cleaned.txt
#^the above was killing screen for some bizarre reason...
cp $SNPS init_snp_lists/${H}.cleaned.txt
###Pass into plink
echo "Preparing for plink"
cut -f 1,2 -d ":" init_snp_lists/${H}.cleaned.txt > init_snp_lists/${H}.cleaned.ids.txt
#
mkdir -p ./pruned_snp_sets/${H}/500kb_0.1r2
 
  plink2 --bfile /work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/ref \
--extract init_snp_lists/${H}.cleaned.ids.txt \
--indep-pairwise 500kb 0.1 \
--out ./pruned_snp_sets/${H}/500kb_0.1r2 \
--exclude /work-zfs/abattle4/ashton/prs_dev/prs_tools/long_range_ld_removal.hg19.tsv


##Extract from GWAS file
mkdir -p gwas_extracts/${H}/500kb_0.1r2
echo "Getting the GWAS stuff..."
python /work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/quickGWASIter.py  --type reg --output gwas_extracts/${H}/500kb_0.1r2/run --gwas_list /work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv --snp_list ./pruned_snp_sets/${H}/500kb_0.1r2.prune.in --extension ".both_sexes.tsv" --gwas_dir /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/unzipped/ #& 

mkdir -p eqtl_extracts/${H}/500kb_0.1r2
#Get the hg38 lookup information- convert hg19 to hg38
awk '(FNR == NR) {arr[$1]; next} ($5 in arr) {print $4":"$2":"$3}' ./pruned_snp_sets/${H}/500kb_0.1r2.prune.in /work-zfs/abattle4/ashton/reference_data/ukbb_1KG_gtex.intersected.no_ambig.no_multi.no_indel.tsv | sed 's/chr//g'  | sort -u > eqtl_extracts/${H}/500kb_0.1r2/500kb_0.1r2.hg38.lookup
#
#
ml python/3.7-anaconda
python3 ${GD}/gtexv8_extract.py eqtl_extracts/${H}/500kb_0.1r2/500kb_0.1r2.hg38.lookup ./eqtl_extracts/${H}/500kb_0.1r2/ --ciseqtls

echo "Data all extracted...."
##Analysis
echo "Run the joint analysis now..."
mkdir -p results/${H}/500kb_0.1r2



#Rscript ${GD}/jointDecomp.R --eqtl_source ./eqtl_extracts/${H}/500kb_0.1r2/ --gwas_source gwas_extracts/${H}/500kb_0.1r2/run.z.tsv --gwas_names gwas_extracts/trait.names.txt --eqtl_names eqtl_extracts/tiss.names.txt --output ./results/${H}/500kb_0.1r2/full_ --na_method 'ZERO'

#Rscript ${GD}/jointDecomp.R --eqtl_source ./eqtl_extracts/${H}/500kb_0.1r2/ --gwas_source gwas_extracts/${H}/500kb_0.1r2/run.z.tsv --gwas_names gwas_extracts/trait.names.txt --eqtl_names eqtl_extracts/tiss.names.txt --output ./results/${H}/500kb_0.1r2/covar_ --na_method 'ZERO' --covar_control

#Rscript ${GD}/jointDecomp.R --eqtl_source ./eqtl_extracts/${H}/500kb_0.1r2/ --gwas_source gwas_extracts/${H}/500kb_0.1r2/run.z.tsv --gwas_names gwas_extracts/trait.names.txt --eqtl_names eqtl_extracts/tiss.names.txt --output ./results/${H}/500kb_0.1r2/maxAvg_ --na_method 'ZERO' --eqtl_selection "MAX_AVG_Z"

Rscript ${GD}/jointDecomp.R --eqtl_source ./eqtl_extracts/${H}/500kb_0.1r2/ --gwas_source gwas_extracts/${H}/500kb_0.1r2/run.z.tsv --gwas_names gwas_extracts/trait.names.txt --eqtl_names eqtl_extracts/tiss.names.txt --output ./results/${H}/500kb_0.1r2/maxAvg_ --na_method 'ZERO' --eqtl_selection "MAX_GENE"

