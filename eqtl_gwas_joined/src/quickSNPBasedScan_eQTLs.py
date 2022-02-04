#11/23
#Goal of this script is to pull out the top snp-gene pair for each gene for each tissue.
#INPUT: Table of snp-gene pairs of pvalues from GTEx (those filtered by YUAN)
#INPUT: list of all SNP-gene pairs we will take under consideration.
#OUTPUT: 3 files with the t-stats, beta, and se of a SNP x Tissue matrix, where entries correspond to the SNP-Gene pair with the strongest signal for a given SNP in a given tissue..

#quick scan
import glob
import sys
lookup_file = "./gene_based_extracts/lookup_snps.hg38.txt"
#TODO: change this to be a pruned set of SNPs input file with the genes. We don't need all these, do we now.:
SNP=1
GENE=0
snp_dict = dict() #dictionary of snps:genes:pvalss
gene_set = dict()
print("Uploading query data...")
with open(lookup_file, 'r') as istream:
  for line in istream:
    dat = line.strip().split()
    if dat[SNP] not in snp_dict:
      snp_dict[dat[SNP]] = dict()
    #snp_dict[dat[0]][dat[1]] = -1 #keep track of all the genes
    gene_set[dat[GENE]] = "" #keep track of all the GENES

yuanfile = "/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/cbset_datasets/input_pairs_fitModel/v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_1_pvalue.txt"
with open(yuanfile, 'r') as istream:
  for line in istream:
    dat = line.strip('\n').split('\t') #just want to strip off the new line, not the whitespace.
    if dat[SNP] in snp_dict:
         snp_dict[dat[SNP]][dat[GENE]] = dat[2:] 

print("Processing results")
ret_set = dict()
n_tiss = 49
count = 0
for snp in snp_dict:
    if "" in snp_dict[snp]:
        continue
    #We now need to iterate through every tissue
    #TODO: run it up to here, check the output
    #For some reason, we keep getting some empty snp-gene pairs, which I don't want.
    ret_set[snp] = dict()
    print(snp)
    for tiss in range(0,n_tiss):
        top_gene = ""
        min_p = 100
        for gene in snp_dict[snp]:
            if (snp_dict[snp][gene][tiss] != "") and (float(snp_dict[snp][gene][tiss]) < min_p):
                top_gene = gene
                min_p = float(snp_dict[snp][gene][tiss])
        if top_gene not in ret_set[snp]:
            ret_set[snp][top_gene] = set()
        ret_set[snp][top_gene].add(tiss)
        print(ret_set[snp])
    count += 1
    if count > 5:
        break
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
