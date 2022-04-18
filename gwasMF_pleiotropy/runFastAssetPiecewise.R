#!/usr/bin/env Rscript
pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr,stats, cowplot, magrittr)
args = commandArgs(trailingOnly=TRUE)
#1 is start SNP
#2 is end SNP

library(ASSET)
start <- as.numeric(args[1])
end <- as.numeric(args[2])
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/quickLoadData.R")
dat.gwas<- quickLoadFactorization("B_SE", "MARCC")
trait.list <- colnames(dat.gwas$X)

ldsc.mat <- as.matrix(fread("/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/ldsc_intercepts.cor.tsv"))
colnames(ldsc.mat) <- gsub(x=colnames(ldsc.mat), pattern = "X", replacement = "")
rownames(ldsc.mat) <- colnames(ldsc.mat)
load("/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/block_0.7.RData")
#run it
nsnps <- (end - start) + 1

p.mat <- matrix(0, nrow = nsnps, ncol = 3) #count 1, count 2, pvalue
for(i in start:end)
{
  tryCatch(
    expr = {
      retsnp <- fast_asset(dat.gwas$vars$ids[i], trait.list, beta.hat = dat.gwas$X[i,], 
                      sigma.hat = 1/dat.gwas$W[i,], Neff = dat.gwas$N[i,2:56], cor = ldsc.mat, block = block)
      
      p.mat[i,1] <- sum(retsnp$Subset.2sided$pheno.1)
      p.mat[i,2] <- sum(retsnp$Subset.2sided$pheno.2)
      p.mat[i,3] <- sum(retsnp$Subset.2sided$pval)
    },
    error = function(e){ 
      message("NA happened :(")
      message(paste0("Error at,", i))
      p.mat[i,1] <- NA
      p.mat[i,2] <- NA
    },
    warning = function(w){
      message("Warning")
      message(paste0("Error at,", i))
      p.mat[i,1] <- NA
      p.mat[i,2] <- NA
    },
    finally = {
      #print(dat.gwas$vars$ids[i])
    }
  )
  
}

o <- data.frame(p.mat) %>% set_colnames(c("cpos", "cneg", "pval")) 
write.table(o, file = paste0("fastAsset_", start, "-", end, ".tsv"), sep = '\t', row.names = FALSE, quote = FALSE)
