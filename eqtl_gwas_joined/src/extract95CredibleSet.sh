#!/bin/bash
#input is list of SNPs in hg37 to query
QUERY=/work-zfs/abattle4/ashton/snp_networks/scratch/eqtl_gwas_joined/fourth_pass_universal/pruned_snp_sets/gwas_centric_1e-5/500kb_0.1r2.prune.in
ODIR="/work-zfs/abattle4/ashton/snp_networks/scratch/eqtl_gwas_joined/fourth_pass_universal/95CredibleSet/gwas1e-5/"
#QUERY=$1
#ODIR=$2
mkdir -p $ODIR
ls -d /work-zfs/abattle4/lab_data/ldsc_qtl_hormozdiari_2018/Annot/*_Analysis_95CredibleSet/ | cut -f 7 -d "/" | sed 's/_95CredibleSet//g' > ${ODIR}/tissue.annotnames.txt
for TISS in `cat ${ODIR}/tissue.annotnames.txt`; do
    echo $TISS
    touch ${ODIR}/${TISS}.annot
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    LOOK="/work-zfs/abattle4/lab_data/ldsc_qtl_hormozdiari_2018/Annot/${TISS}_95CredibleSet/${TISS}_95CredibleSet.${i}.annot"
    echo $LOOK
    cat ${ODIR}/${TISS}.annot <(awk '(FNR == NR) {arr[$1];next} ($1":"$2 in arr) {print $0}' $QUERY $LOOK | sort -h -k2) > t && mv t ${ODIR}/${TISS}.annot
    done
done

