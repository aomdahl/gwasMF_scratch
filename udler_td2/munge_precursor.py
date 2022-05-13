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
snp_options = ["MARKER"] 
a1_options = ["RISK_ALLELE"] #effect allele
a2_options = ["OTHER_ALLELE"] #alternate allele
sum_stats = ["MAINEFFECTS"]
se_options = ["MAINSE"]
pval_options = ["MAINP"]
labels = {"n" : n_options, "snp" : snp_options, "a1" : a1_options, "a2" :a2_options, "ss" : sum_stats, "se" :se_options, "p" : pval_options}
import sys
import argparse
import gzip
import os
from subprocess import check_call
import subprocess
import shlex
import csv
import numpy as np

which = lambda lst:list(np.where(lst)[0])


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
        return gzip.open(path, 'rb')
    else:
        return open(path, 'r')

def outFileContext(path):
    if path[-3:] == ".gz" or path[-4:] == ".bgz":
        return gzip.open(path, 'wt')
    else:
        return open(path, 'w')
    
def detectDelimiter(line):
    #from https://stackoverflow.com/questions/3952132/how-do-you-dynamically-identify-unknown-delimiters-in-a-data-file
    sniffer = csv.Sniffer()
    return sniffer.sniff(line).delimiter

#File handler- if gz or not.
def fh(st, p):
    try:
        if ".bgz" in  p or ".gz" in p:
            return st.decode(errors='replace')
        else:
            return st
    except TypeError:
        print(p)
        print(str(p))
        input()
        return p

def labelSpecify(header, checkfor, labels):
    """
    Specifies which column maps to which thing in cases when munge won't pick it up
    @header- the LIST of header items
    @checkfor: which item we are looking for
    @labels: the map of labels to lookup options.
    """
    str_map = {
        "snp" : " --snp ",
        "n" : " --N-col ",
        "a1": " --a1 ",
        "a2": " --a2 ",
        "ss": " --signed-sumstats "
    }
    relevant_col = which([x.upper() in labels[checkfor] for x in header])
    if(len(relevant_col) > 0):
        ret_snp = str_map[checkfor] + str(header_dat[snp_col[0]]) + " "
    else: 
        ret_snp = ""
    return ret_snp

def filePeek(readin):
    """
    This looks at the summary stats file and determines if munge_stats will be compatible with it.
    If not, it makes necessary changes or notifies the user. Yup.
    @return cleanup_protocol: instructions for cleaning up the file, if any.
    @return ret_n the string to append to the munge command to deal with sample sizes
    """
    fpath = readin[0]
    pheno = readin[1]
    
    ret_n = ""
    try:
        with openFileContext(fpath) as istream:
            header = fh(istream.readline(), fpath)
            header = header.strip()
            first = fh(istream.readline(),fpath)
            first = first.strip().split()
            delim  = detectDelimiter(header)

            if delim == ",":
                cleanup_protocol = "CSV"
            
            if ("PaxHeader" in header) or ("GIANT" in fpath):
                #This is a GIANT file, need to do cleanup
                cleanup_protocol = cleanup_protocol + "GIANT"
                return cleanup_protocol, ret_n

            header_dat = header.split()
            if header_dat[0] == "variant" and header_dat[-1] == "pval": #we suspect this is UKBB full format file.
                if isVariantID(first[0]):
                    print("The current file at", fpath, "appears to be a full Neale Lab UKBB file. We don't recommend using this at this stage, use the LDSC pre-formatted one available onlinel")
                    cleanup_protocol = cleanup_protocol + "UKBB"
                    return cleanup_protocol, ret_n
                    #need to deal with UKBB cases of both continuous and case-control

            if header_dat[0].upper() == "CHR" and header_dat[1].upper() == "POS":
                cleanup_protocol = cleanup_protocol + "TO_RSID"
            if ":" in header_dat[0]:
                cleanup_protocol = cleanup_protocol + "TO_RSID_1"

            #Get the sample size if its there...
            ret_n =""
            if any([x.upper() in labels["n"] for x in header_dat]):
                next #n will be detected by munge
            else:
                #N will not be detected by it...
                if any(["N_SAMPLES" in x.upper() for x in header_dat]):
                #Specify the column for N
                    ret_n = labelSpecify(header_dat, "n", labels)
                else: #specify the input from the file,
                    ret_n = " --N " + dat[3]

            #Check for the marker data
                ret_a1 = labelSpecify(header_dat, "a1", labels)
                ret_a2 = labelSpecify(header_dat, "a2", labels)
                ret_ss = labelSpecify(header_dat, "ss", labels)
                ret_snp = labelSpecify(header_dat, "snp", labels)

    except FileNotFoundError:
        print("Specified path does not exist")
        cleanup_protocol = "ERROR"
    
    if cleanup_protocol == "": #nothing to change.
        cleanup_protocol = "NONE"
    return cleanup_protocol, [ret_a1, ret_a2, ret_ss, ret_snp, ret_n] 

def makeInterFileName(fpath, protocol):
    """
    If no change to be made, don't modify the file.
    """
    #os.path.basename()
    if "TO_RSID" in protocol or "GIANT" in protocol :
        return os.path.splitext(fpath)[0] + ".INTERMEDIATE" + os.path.splitext(fpath)[1]
    else:
        return fpath

def buildOutName(outpath, inpath, inpheno):
    """
    Specifies the name of the output file based on the output directory, the phenotype name, etc.
    @param outpath: the provided output directory
    @param inpath: the input file info
    """
    fn = os.path.basename(inpath)
    with_ext = os.path.splitext(fn)

    return outpath + with_ext[0] + "." + inpheno

def doCleanup(readin, protocol, rsids):
    """
    Performs the cleanup step on files that need it
    UKBB, we shouldn't do anything
    Giant: remove teh first few lines, and get the header, keep the rest of the file the same
    TO_RSID: convert the chr/things to rsids; make sure to check using the right genome build (!)
    """
    fpath = readin[0]
    intermediate_file_name = makeInterFileName(fpath, protocol)
    delim = ""
    if "GIANT" in protocol :
        #remove the first few lines of the file
        #done_header="\"MarkerName\\tChr\\tPos\\tAllele1\\tAllele2\\tFreqAllele1HapMapCEU\\tb\\tse\\tp\\tN\""
        done_header="MarkerName\tChr\tPos\tAllele1\tAllele2\tFreqAllele1HapMapCEU\tb\tse\tp\tN\n"
        print_switch = False
        with openFileContext(fpath) as istream:
            with outFileContext(intermediate_file_name) as ostream:
                for i, line in enumerate(istream):
                    line = fh(line.strip(), fpath)
                    
                    if print_switch:
                        ostream.write(line + '\n')
                    if "MarkerName" in line:
                        ostream.write(done_header)
                        print_switch = True
           
    if "TO_RSID" in protocol:
        #If its TO_RSID, its a 2 column thing. if its TO_RSID_1 its a 1 column thing.
        dropped_snps = 0
        T='\t'
        print("Converting SNP IDs to RSIDs, constraining to only Hapmap3 non-HLA SNPs")
        with openFileContext(fpath) as istream:
            with outFileContext(intermediate_file_name) as ostream:
                for i, line in enumerate(istream):
                    line = fh(line.strip(), fpath)                     
                    if i == 0: 
                        delim = detectDelimiter(line)
                        ostream.write("SNP" + T + T.join(line.split(delim)[1:]) + '\n')

                    else:
                        line = line.split(delim)
                        if protocol == "TO_RSID":
                            address = (line[0] + ":" + line[1]).replace("chr", "")
                        else:
                            address = line[0].replace("chr", "")
                        if address in rsids:
                            snp = rsids[address]
                            ostream.write(snp + T + T.join(line[1:]) +  '\n')
                        else:
                            dropped_snps += 1
                            continue
        print(dropped_snps, "SNPs omitted.")
    #For all remaining protocols, we don't do anything.                    
        #determine which columnes have the info we need, look it up, and convert it.
    
    if "CSV" in protocol:
        if "TO_RSID" in protocol: #we have had to modify other things too..
            #nothing to change, it was edited above
            break
        else:
            print("Changing CSV to TSV...")
            with openFileContext(fpath) as istream:
                with outFileContext(intermediate_file_name) as ostream:
                    for i, line in enumerate(istream):
                        line = fh(line.strip(), fpath)
                        ostream.write(T.join(line.split()) + '\n')                    
            input(("Check comma conversion"))
        

    return intermediate_file_name


def importRSDict(p):
    """
    By default uses a list of no hLA region hm3 snps
    """
    ret_dict = dict()
    with open(p, 'r') as istream:
        for line in istream:
            line = line.strip().split()
            #ret_dict[line[0] + ":" line[1]] = line[2]
            ret_dict[line[0]] = line[1]
    return ret_dict


parser = argparse.ArgumentParser(description = "Quick script to get many GWAS summary stats converted into an LDSC-friendly format. A nice wrapper for ldsc's munge script.")
parser.add_argument("--study_list", help = "list of studies to work on. Columns are study path, phenotype name, genome build, N samples. Last 2 columns not required, in which case we assume hg37 and that N is in the data. Note that if N is in the data, that is prioritized above a specified N.")
parser.add_argument("--rsid_ref", help = "path to the RSID reference")
parser.add_argument("--ldsc_path", help = "path to LDSC code", default="/work-zfs/abattle4/ashton/genomics_course_2020/project_2/ldsc/")
parser.add_argument("--output", help = "output path")
args = parser.parse_args()

#get the hg37 build dict in memory....
build_list = {"hg37":importRSDict("/work-zfs/abattle4/ashton/reference_data/hm3_nohla.snpids"), "hg38":None, "hg36": None}
with open(args.study_list ,'r') as istream:
    with open(args.output + "munge_sumstats.all.sh", 'w') as ostream:
        for line in istream:
            dat = line.strip().split()
            print("Currently processing file ", dat[0], "(", dat[1], ")")
            #Check the first lines of the file, see if its UKBB template, if it has N, etc.
            #Check to see if its got RSIDs or now...
            fix_instructions, fixed_labels = filePeek(dat)
            #detect the genome build if needed....
            if (len(dat) > 2) and (dat[2] != "" )and ("TO_RSID" in fix_instructions):
                print("Input specifies", dat[2], "build.")
                if(dat[2] != "hg37" and dat[2] != "hg19"):
                    print("not yet implemented")
                    continue
                else:
                    build_dict = build_list["hg37"]
            else:
                print("Assuming hg37 build, using hapmap3 SNPS only...")
                build_dict = build_list["hg37"]
            intermediate_file = doCleanup(dat, fix_instructions,build_dict)
            #build commands for it
            if fix_instructions != "UKBB" and fix_instructions != "ERROR":
                lc = "python2 " + args.ldsc_path + "munge_sumstats.py --sumstats " + intermediate_file + " --out " + buildOutName(args.output, dat[0], dat[1]) + "".join(fixed_labels)
                ostream.write(lc + '\n')
            print()
        

       



"""
                        #this is the one we want
                        command = "zcat " + fpath +  " | awk '(FNR > " + str(i + 1) + ") {print $0}'"
                        #command = "cat <(echo -e " + done_header + ") <(zcat /work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/data/GIANT_2015_WHR_COMBINED_EUR.txt.gz | awk '(FNR > " + str(i + 1) + ") {print $0}')"
                        #issue- subprocess doesn't like piping :(
                        # | gzip > " + intermediate_file_name 
                        print(intermediate_file_name)  
                        print(shlex.split(command))
                        output = gzip.open(intermediate_file_name, 'wt')
                        output.write(done_header)
                        subprocess.Popen(shlex.split("zcat " + fpath), stdout = subprocess.PIPE)
                        subprocess.Popen()
                        #subprocess.Popen(shlex.split(command), stdout = output)
                        output.close
                        #subprocess.Popen()
                        ml python/3.7-anaconda
                        python munge_precursor.py --study_list pheno_try.txt --output ./ldsr_all/
"""




