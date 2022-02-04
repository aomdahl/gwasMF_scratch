#!/usr/bin/env python3
import glob
import sys
import os
#Quickly extract snps from GTEX v8 
#@arg dir_path: directory

def writeFile(fl, args):
    handle = args[2] + os.path.basename(fl)[:-4] + ".extract.tsv"
    print(handle)
    with open(fl, 'r') as istream:
        with open(handle, 'w') as ostream:
            counter = 0
            for l in istream:
                counter += 1
                line = l.strip().split()
                snp = line[1].replace("_", ":").replace("chr", "")[:-4] #this will give 1:1231:G:A
                if snp in lookup:
                    #print it right away, avoid storing then printing
                    #columns of interest are all of them
                    ostring = line[0] + "\t" + snp + "\t" + "\t".join(line[2:]) + "\n"
                    #print(ostring)
                    ostream.write(ostring)
                if(counter > 10000000):
                    print("ten mil lines in")
                    counter = 0











if sys.argv[1] == "-h" or sys.argv[1] == "--help":
    print("this utility will quickly extract a set of snps from all gtex v8 tissue files. Output will be given in order it appears in GTEx files")
    print("Arguments are \n 1) the list of SNP ids to extract 2) output prefix 3) Which source (spQTLs or ciseQTLs) \n 4) which tissue to start from (optional)")
    sys.exit()
lookup = dict()
#read in list of SNPs to include
with open(sys.argv[1], 'r') as istream:
    for l in istream:
        lookup[l.strip()] =""
if sys.argv[3] == '-splicing':
    print("doing the splicing thing....")
    files = glob.glob("/work-zfs/abattle4/lab_data/GTEx_v8/sqtl/GTEx_Analysis_v8_sQTL_all_associations/*allpairs.txt")
else:
    print('Expression it is')
    files = glob.glob("/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/*allpairs.txt")
#allow to start farther along- at another file, if it got interrupted.
if len(sys.argv) == 4 :
    start_at = 0
else:
    start_at = [idx for idx, s in enumerate(files) if sys.argv[4] in s][0]
    print("Starting at " + str(start_at))
    print(files[start_at])


#header style
'''
gene_id variant_id  tss_distance    ma_samples  ma_count    maf pval_nominal    slope   slope_se
'''
from joblib import Parallel, delayed
import multiprocessing
query_list = files[start_at:]
num_cores = multiprocessing.cpu_count()
if num_cores > 20:
    num_cores = 20
Parallel(n_jobs = num_cores)(delayed(writeFile)(f,sys.argv) for f in query_list)
"""

for fl in files[start_at:]:
    print(fl)
    handle = sys.argv[2] + os.path.basename(fl)[:-4] + ".extract.tsv"
    print(handle)
    with open(fl, 'r') as istream:
        with open(handle, 'w') as ostream:
            counter = 0
            for l in istream:
                counter += 1
                line = l.strip().split()
                snp = line[1].replace("_", ":").replace("chr", "")[:-4] #this will give 1:1231:G:A
                if snp in lookup:
                    #print it right away, avoid storing then printing
                    #columns of interest are all of them
                    ostring = line[0] + "\t" + snp + "\t" + "\t".join(line[2:]) + "\n"
                    #print(ostring)
                    ostream.write(ostring)
                if(counter > 10000000):
                    print("ten mil lines in")
                    counter = 0

"""
