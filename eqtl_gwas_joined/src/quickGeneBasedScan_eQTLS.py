#read in the list, copy from elsewhere
#INPUT: Table of snp-gene pairs of pvalues from GTEx (those filtered by YUAN)
#INPUT: list of all SNP-gene pairs we will take under consideration.
#OUTPUT: 3 files with the t-stats, beta, and se of a Gene x Tissue matrix, where entries correspond to the the Gene x SNP pair with the lowest p-value for that Gene and Tissue.


#quick scan
import glob
import sys
lookup_file = "./gene_based_extracts/lookup_snps.hg38.txt"
#ofile = "quick_test.txt"
#ofile = sys.argv[2]
gene_dict = dict() #dictionary of genes:snps:pvals
snp_set = dict()
print("Uploading query data...")
with open(lookup_file, 'r') as istream:
  for line in istream:
    dat = line.strip().split()
    if dat[0] not in gene_dict:
      gene_dict[dat[0]] = dict()
    #gene_dict[dat[0]][dat[1]] = -1 #keep track of all the genes
    snp_set[dat[1]] = "" #keep track of all the SNPs
#print(gene_dict)
#okay, thinking this through the SNPs don't actually matter, its only the GENEs that I care about here.
#unless we want to limit it to SNPs that COULD BE in both, which maybe is a good idea....?
#not sure that it really matters. We are limiting the UKBB set though, so might be best to be fair?
#doesn't need to be a perfect analysis, just run it and see what you get.
GENE=0
SNP=1
yuanfile = "/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_1_pvalue.txt"
with open(yuanfile, 'r') as istream:
  for line in istream:
    dat = line.strip('\n').split('\t') #just want to strip off the new line, not the whitespace.
    if dat[GENE] in gene_dict:
         gene_dict[dat[GENE]][dat[SNP]] = dat[2:] 
        #if dat[SNP] in gene_dict[dat[GENE]]:
        #    #gene_dict[dat[GENE]] = dict()
        #    #print("Euston, we have a problem")
        #    gene_dict[dat[GENE]][dat[SNP]] = dat[2:] #all the rest.
            

print("Processing results")
ret_set = dict()
n_tiss = 49
for gene in gene_dict:
    #We now need to iterate through every tissue
    ret_set[gene] = dict()
    print(gene)
    for tiss in range(0,n_tiss):
        top_snp = ""
        min_p = 100
        for snp in gene_dict[gene]:
            if (gene_dict[gene][snp][tiss] != "") and (float(gene_dict[gene][snp][tiss]) < min_p):
                top_snp = snp
                min_p = float(gene_dict[gene][snp][tiss])
        #if top_snp == "":
        #    print("No matches found?")
        #    print(gene_dict[gene])
        #    input()
        if top_snp not in ret_set[gene]:
            ret_set[gene][top_snp] = set()
        ret_set[gene][top_snp].add(tiss)
        #print(tiss, gene, top_snp, min_p)
        #input()
#what remains to be done:

ifile="/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_1_se.txt"
#specify default value, 1 or 0
def readYuanMat(ifile, def_val):
    print_dict = dict()
    pos_dict = dict()
    with open(ifile, 'r') as istream:
        for line in istream:
            dat = line.strip('\n').split('\t')
            if (dat[GENE] in ret_set) and (dat[SNP] in ret_set[dat[GENE]]):
                if dat[GENE] not in print_dict:
                    print_dict[dat[GENE]] = [def_val] * 49
                    pos_dict[dat[GENE]] = [""] * 49
                for pos in ret_set[dat[GENE]][dat[SNP]]:
                    #print(dat)
                    print_dict[dat[GENE]][pos] = dat[pos + 2]
                    pos_dict[dat[GENE]][pos] = dat[SNP]
                    #print("Inserting at", dat[SNP],dat[GENE])
                    #print("Tissue", pos)
                    #print("Value", dat[pos + 2])
                    #input()
    return print_dict, pos_dict

se, se_chr = readYuanMat(ifile, 1)
slope, slope_chr = readYuanMat("/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_1_slope.txt", 0)
import pandas as pd
se_pd = pd.DataFrame.from_dict(se) #,orient = 'index')
se_pd = se_pd.apply(pd.to_numeric)
beta_pd = pd.DataFrame.from_dict(slope) #,orient = 'index')
beta_pd = beta_pd.apply(pd.to_numeric)
z_pd = beta_pd/se_pd

#get the addresses:
se_add = pd.DataFrame.from_dict(se_chr)
beta_add = pd.DataFrame.from_dict(slope_chr)

z_pd.to_csv("./gene_based_extracts/eqtl_z_scores.csv",index=False)
se_add.to_csv("./gene_based_extracts/se_snps.csv", index=False)
beta_add.to_csv("./gene_based_extracts/beta_snps.csv", index=False)
