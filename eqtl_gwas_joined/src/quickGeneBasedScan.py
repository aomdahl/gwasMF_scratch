#quick scan
import glob
import sys
lookup_file = "./gene_based_extracts/lookup_snps.txt"
#ofile = "quick_test.txt"
ofile = sys.argv[2]
lookup_dict = dict()
lookup_set = dict()
print("Uploading query data...")
with open(lookup_file, 'r') as istream:
  for line in istream:
    dat = line.strip().split()
    if dat[0] not in lookup_dict:
      lookup_dict[dat[0]] = list()
    lookup_dict[dat[0]].append(dat[1]) #keep track of all the genes
    lookup_set[dat[1]] = "" #keep track of all the SNPs

print("Parsing summary stats...")
#this step is surprisingly fast
#gwas_file="/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/unzipped/21001_irnt.gwas.imputed_v3.both_sexes.tsv"
gwas_file = glob.glob("/work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/unzipped/" + sys.argv[1] + "*.both_sexes.tsv")[0]
print(gwas_file) 
#gwas_file = sys.argv[1]
gwas_dict = dict()
count = 1
with open(gwas_file, 'r') as istream:
  for line in istream:
    dat = line.strip().split()
    STAT=9
    if len(dat) == 12:
      STAT = 10
    
    if dat[0] in lookup_set:
      gwas_dict[dat[0]] = float(dat[STAT]) #get the zscore
    count += 1
    #if count % 1000000 == 0:
    #  print(count)
#print("Done")

print("Parsing by gene...")

with open(ofile, 'w') as ostream:
  for gene in lookup_dict:
    #minp = 100.0
    #min_snp = ""
    max_t = 0.0
    max_snp = ""
    for snp in lookup_dict[gene]:
      #print("Current SNP", snp)
      #print("Current z score", gwas_dict[snp])
      #if float(gwas_dict[snp]) < minp:
      if abs(gwas_dict[snp]) > abs(max_t):
        #minp = float(gwas_dict[snp])
        #min_snp = snp
        max_snp = snp
        max_t = gwas_dict[snp]
    #print(gene, min_snp, minp)
    ostream.write(gene + "\t"  + max_snp + "\t" +  str(max_t) + "\n")
