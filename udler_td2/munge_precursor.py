""" 
for f in files
    read in first few lines
        identify headers, and each assignment
        if is chr:00 and hg37:
            read in, do conversion lookup
        if effect/alt is lowercase, make uppercase
        call to mung 
"""
n_options =['N','NCASE','CASES_N','N_CASE', 'N_CASES','N_CONTROLS','N_CAS','N_CON','N_CASE','NCONTROL','CONTROLS_N','N_CONTROL']

import sys
import argparse
import numpy as np
import glob
import pickle
import gzip


which = lambda lst:list(np.where(lst)[0])

def buildOutName(infoline):
    """
    Specifies the name of the output file based on the output directory, the phenotype name, etc.
    
    """
    return ""
def isVariantID(text):
    """
    Helper function to detect if this has addresses, not rsids, which we want.
    """
    split = text.split(":")
    if (split[0] in [str(x) for x in list(range(0,24))]) or ("CHR" in split[0].upper()): 
        return True
    else:
        return False


def openFileContext(path):
    if path[-3:] == ".gz" or path[-4:] == ".bgz":
        print("A bzip file...")
        return gzip.open(path, 'rb')
    else:
        return open(path, 'r')

def filePeek(readin):
    """
    This looks at the summary stats file and determines if munge_stats will be compatible with it.
    If not, it makes necessary changes or notifies the user. Yup.
    @return cleanup_protocol: instructions for cleaning up the file, if any.
    @return ret_n the string to append to the munge command to deal with sample sizes
    """
    fpath = readin[0]
    pheno = readin[1]
    cleanup_protocol = "NONE"
    ret_n = ""
    with openFileContext(fpath) as istream:
        header = istream.readline().decode(errors='replace')
        header = header.strip()
        first = istream.readline().decode(errors='replace')
        first = first.strip().split()

        
        if ("PaxHeader" in header) or ("GIANT" in fpath):
            #This is a GIANT file, need to do cleanup
            cleanup_protocol = "GIANT"
            return cleanup_protocol, ret_n

        header_dat = header.split()
        print(header_dat)

        if header_dat[0] == "variant" and header_dat[-1] == "pval": #we suspect this is UKBB full format file.
            if isVariantID(first[0]):
                print("The current file at", fpath, "appears to be a full Neale Lab UKBB file. We don't recommend using this at this stage, use the LDSC pre-formatted one available onlinel")
                cleanup_protocol = "UKBB"
                return cleanup_protocol, ret_n
                 #need to deal with UKBB cases of both continuous and case-control

        if header_dat[0].upper() == "CHR" and header_dat[1].upper() == "POS":
            cleanup_protocol = "TO_RSID"
        #Get the sample size if its there...
        ret_n =""
        if any([x.upper() in n_options for x in header_dat]):
            next #n will be detected by munge
        else:
            #N will not be detected by it...
            if any(["N_SAMPLES" in x.upper() for x in header_dat]):
               #Specify the column for N
                n_col = which(["N_SAMPLES" in x.upper() for x in header_dat])
                ret_n ="  --N-col " + header_dat[n_col]
            else: #specify the input from the file,
                ret_n = " --N " + dat[2]
    return cleanup_protocol, ret_n

def doCleanup(readin, protocol):
    fpath = readin[0]
    if protocol == "UKBB":
        return
    if protocol == "GIANT":
        #remove the first few lines of the file
        done_header="MarkerName	Chr	Pos	Allele1	Allele2	FreqAllele1HapMapCEU	b	se	p	N"
        with openFileContext(fpath) as istream:
            for i, line in enumerate(istream):
                print(i)
                line = istream.readline().decode(errors='replace')
                if "MarkerName" in line:
                    #this is the one we want
                    run = "cat <(echo -e " + done_header + ") <(zcat " + fpath +  " | tail -n +" + str(i + 1) + ") | gzip > outfile.txt"
                    print(run)
                    return
            
    if protocol == "TO_RSID":
        return
        #determine which columnes have the info we need, look it up, and convert it.


parser = argparse.ArgumentParser(description = "Quick script to get many GWAS summary stats converted into an LDSC-friendly format. A nice wrapper for ldsc's munge script.")
parser.add_argument("--study_list", help = "list of studies to work on. Columns are study path, phenotype name, N samples. Last column not required.")
parser.add_argument("--rsid_ref", help = "path to the RSID reference")
parser.add_argument("--ldsc_path", help = "path to LDSC code", default="/work-zfs/abattle4/ashton/genomics_course_2020/project_2/ldsc/")
parser.add_argument("--pickle", default = "", help = "Path to a pickle file.")
parser.add_argument("--output", help = "output path")
parser.add_argument("--gwas_dir", help = "gwas directory to look in.")
parser.add_argument("--type", help = "if its ldsc or not", default = "NOT")
parser.add_argument("--extension", help = "File extension:q to grab")
args = parser.parse_args()

with open(args.study_list ,'r') as istream:
    for line in istream:
        dat = line.strip().split()
        #Check the first lines of the file, see if its UKBB template, if it has N, etc.
        #Check to see if its got RSIDs or now...
        fix_instructions, ncount = filePeek(dat)
        #build commands for it
        print(args.ldsc_path + "munge_sumstats.py --sumstats " + dat[0] + " --out " + args.output+ buildOutName(dat) + ncount)
        print(fix_instructions)
        doCleanup(dat, fix_instructions)

       



