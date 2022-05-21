#!/bin/bash
#Workflow to run LDSC on everything

ml python/3.7-anaconda
cd /work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2

#prep the summary stats for LDSC analysis.
python src/munge_precursor.py --study_list pheno_try.sansUKBB.txt --merge_alleles /work-zfs/abattle4/ashton/reference_data/hm3_snps_ldsc_ukbb.tsv --output /work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/may16
bash may16munge_sumstats.all.sh 

#Generate the commands to run LDSC pairwise on each set of traits
bash src/all_pairwise_r2g.sh ldsr_all/pairwise.traits.txt /work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/ldsr_results/ 10

#Drop the jobs on slurm
for i in {0..9}; do
cat /work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/run_scripts/slurm_template.sh /work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/ldsr_results/_ldsc.run.${i}.sh > run_scripts/slurm.ldsc.run.${i}.sh
sbatch run_scripts/slurm.ldsc.run.${i}.sh
done

#once the jobs are completed, extract the relevant information.
P=/work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/src
for i in ldsr_results/*.log; do
    echo $i
    python $P/parseLDSCr2g.py --input_file $i --output ldsr_results/
    b=`fileName $i`
    grep -A 1000 "gcov_int_se" $i  | grep -v "Analysis finished" | grep -v "Total time" > ldsr_results/$b.rg_report.tsv
done
