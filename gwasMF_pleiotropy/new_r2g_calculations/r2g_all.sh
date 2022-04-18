#!/bin/bash
set -e
ODIR=/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/new_r2g_calculations/
FILE=runfiles.txt
#Get of list of files to run on...
N=`wc -l $FILE| awk '{print $1}'`
rm runames.tmp
echo $N
for ((i=1;i<=$N;i++)); do
     #drop the current one we want
     sed -e ${i}d $FILE | tr '\n' ',' | sed 's/,$//g' >> runames.tmp
     echo "" >> runames.tmp
done

#Get the output file names based on the reference GWAS
rm runids.tmp
while read p; do
     basename $p | cut -f 1 -d "." >> runids.tmp
done < $FILE



for ((i=1;i<=$N;i++)); do
  
  QUERY=`sed -n ${i}p runames.tmp`
  P=/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/
    ID=`sed -n ${i}p runids.tmp`
cd /work-zfs/abattle4/ashton/reference_data/ldsc_ref/
   python  /work-zfs/abattle4/ashton/genomics_course_2020/project_2/ldsc/ldsc.py \
    --rg $QUERY \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out ${ODIR}/${ID}
    cd -
done < runames.txt
