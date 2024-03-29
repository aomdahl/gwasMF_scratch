#!/bin/bash
#To run a simulation, create your desired yaml file and specify an output directory

ml anaconda
    conda activate renv
    #Just for the first case, create 10 replicates of X, SE to run from.
    YAML=$1
    #YAML=simulating_factors/udler_based_500/udler_1_3block-corr_realistic-high-overlap/udler_1_3block-corr_realistic-high-overlap.yml 
    ODIR=$2
    #ODIR=simulating_factors/udler_based_500/udler_1_3block-corr_realistic-high-overlap/
    NITER=`grep "iter" $YAML | cut -f 2 -d ","` #given in the yaml file

        #Now actually run the factorization.
    mkdir -p ${ODIR}/factorization_results
    MNAMES=`grep "test_methods" $YAML | cut -f 2 -d "," | sed 's/:/,/g'`
    K=`grep "K" $YAML | cut -f 2 -d ","`
    for i in $(eval echo "{1..$NITER}"); do
        Rscript /home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/matrix_factorization.R --se_data ${ODIR}/sim${i}.std_error.txt --beta_data ${ODIR}/sim${i}.effect_sizes.txt --seed ${i} --outdir ${ODIR}/factorization_results/sim${i}.whitened. --only_run $MNAMES --K $K --no_plots --C ${ODIR}/sim${i}.c_matrix.txt
    done

#assess the performance of this...
Rscript src/evaluateSims.R --output ${ODIR}/factorization_results/summary --plot --yaml  $YAML --sim_path ${ODIR}/factorization_results/ -w

