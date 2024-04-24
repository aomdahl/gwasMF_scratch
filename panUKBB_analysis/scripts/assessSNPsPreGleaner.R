#!/usr/bin/env Rscript
pacman::p_load(magrittr, dplyr, data.table)

setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization")
ref.lookup <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz") 
p.vals <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.p.tsv")

#Check the 1e-5 threshold
meet.thresh <- apply(p.vals[,-1],1 , function(x) sum(as.numeric(x) < 1e-5,na.rm = TRUE))
bad.snps <- p.vals[which(meet.thresh < 1), 1] %>% print()
if(all(is.na(p.vals[which(meet.thresh < 1),-1])))
{
  message("All SNPs PASS THE P-VALUE threshold")
}


bad.maf <- ref.lookup %>% filter(rsid %in% p.vals$ids) %>% filter(MAF < 0.01)
if(nrow(bad.maf)==0)
{
  message("All SNPs pass MAF threshold")
}



#Missing by study
missing.per.study <- apply(p.vals[,-1],2,function(x) sum(is.na(x)))
summary(missing.per.study)/nrow(p.vals)

if(any(missing.per.study >  nrow(p.vals[,-1])/3))
{
  message("A total of ", sum(missing.per.study >  (nrow(p.vals[,-1])/3)), " studies are missing data for at least 1/3 of SNPs. This should be dropped")
}


#Missing by SNP
missing.per.snp <- apply(p.vals[,-1],1,function(x) sum(is.na(x)))
summary(missing.per.snp)/ncol(p.vals[,-1])
if(any(missing.per.snp >  ncol(p.vals[,-1])/2))
{
  message("A total of ", sum(missing.per.snp >  (ncol(p.vals[,-1])/2)), " SNPs are missing data for over half of studies. These should be dropped")
}

#All unique?

if(length(p.vals$ids) == length(unique(p.vals$ids)))
{
  message("All SNPs are unique")
}

#write out an image
par(mfrow=c(2,2))
#Plot the distribution of # signifcant
hist(meet.thresh,xlab="# SNPs per study at p < 1e-5", main = "Number of sig. SNPs per study", breaks = 30)
#MAF distribution
hist((ref.lookup %>% filter(rsid %in% p.vals$ids))$MAF, xlab = "MAF", main = "Most SNPs have lower MAF")
#Missing per SNP
hist(log10(sort(missing.per.snp)/ncol(p.vals[,-1])), 
     xlab = "Log10-proportion of studies missing for a SNP", main="Proportion missing studies per SNP (log10 scale)", breaks = 50);abline(v=log10(0.5), col="red")
text(c(log10(0.5)),c(200), c("50% missing"),
     cex = 0.88, pos = 2, col = "red")

#missing per study
hist(missing.per.study/nrow(p.vals), breaks = 30, main="Missing SNPs per study", xlab = "Proportion of SNPs missing per study")

