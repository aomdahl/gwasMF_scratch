args = commandArgs(TRUE)
print(args)
set.seed(args[1])
seedi = args[1]
herit_a=0.1
herit_b=0.05
message("Seed set to ", args[1])
message("this version of the script simply saves out all the GWAs summary stats, doesn't actually do any analysis.")
message("herit a: ", herit_a)
message("herit b: ", herit_b)

full_ret=FALSE
if(args[2] == "1")
{
  message("Getting the full report")
  full_ret=TRUE
}
message("Number of simulations to run:", args[3])
n.iter <- as.numeric(args[3])
pacman::p_load(stats, data.table, magrittr, dplyr, tidyr, ggplot2, stringr, fossil)
source("/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/helper_functions_plieo.R")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
#load("/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/.RData")

#Generate genotypes for 3 populations...
## These are based on MAFs from the thousand genomes, subsetted to a pruned set used in real analysis
thou.g.mafs <- fread("/scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/plink.frq")
#thou.g.pruned <- fread("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/500kb.0.04r2.prune.in", header = F)
thou.g.pruned <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark/ukbb_benchmark.250kb.0.2r2.prune.in", header = F)

thou.g.mafs.sub <- thou.g.mafs %>% filter(SNP %in% thou.g.pruned$V1)
thou.g.mafs.use <- thou.g.mafs.sub
rm(thou.g.mafs)

#First- create 3 LARGE cohorts of 50,000 individuals each.

N_tot <- 45000
#1
genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop1.geno <- t(do.call("cbind", genotypes.all))
#2
genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop2.geno <- t(do.call("cbind", genotypes.all))

#3
genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop3.geno <- t(do.call("cbind", genotypes.all))
rm(genotypes.all)


#Simulation settings
rep.dat <- list()
message("Starting sim iters...")
n.snps <- 1000
snp.weights <- list()
snp.weights[[1]] <- c(rnorm(200,sd=sqrt(herit_a/200)), rep(0,800))
snp.weights[[2]] <- c(rep(0,200), rnorm(25,sd=sqrt(herit_b/25)),rep(0,775))
snp.weights[[3]] <- c(rep(0,10), rep(0,900),rep(0.07,90))
stopifnot(length(snp.weights[[3]]) == n.snps)
stopifnot(length(snp.weights[[2]]) == n.snps)
stopifnot(length(snp.weights[[1]]) == n.snps)
for(rep in 1:n.iter)
{
  start = Sys.time() 
 message("on iteration ", rep, " of ", n.iter)
  #Simulate new phenotypes
  overlap.list <- c(0,1000,5000,7500,10000,12500,15000)
  pop1.null <- list()#genotypes.tab
  y1 <- phenotypeBuilder(snp.weights, pop1.geno[,1:n.snps], 1, type = "1+2=3")
  y2 <- phenotypeBuilder(snp.weights, pop2.geno[,1:n.snps], 1, type = "1+2=3")
  y3 <- phenotypeBuilder(snp.weights, pop3.geno[,1:n.snps], 1, type = "1+2=3")
  pop2.null <- list() #genotypes.tab.2
  pop3.null <- list()
  
  n.fixed = 15000
  #start with 2 groups for simplicity in overlapping
  #TODO- report the number of overlapping samples (max = n.fixed, overlap written out each time)
  #TODO- report the correlation in the overlapping phenotype samples
  for(i in overlap.list)
  {
    print(i)
    
    first <- relatedOverlapBETTER(n.snps, i, 3, y1, N.tot = N_tot, 
                                           genotypes=pop1.geno[,1:n.snps], fixed.count = n.fixed,ret.all.stats=full_ret)

   second <- relatedOverlapBETTER(n.snps, i, 3, y2, N.tot = N_tot, 
                                           genotypes=pop2.geno[,1:n.snps], fixed.count = n.fixed,ret.all.stats=full_ret)
    
    third <- relatedOverlapBETTER(n.snps, i, 3, y3, N.tot = N_tot, 
                                           genotypes=pop3.geno[,1:n.snps], fixed.count = n.fixed,ret.all.stats=full_ret)
    if(!full_ret)
    {
      write.table(data.frame(round(cbind(first, second, third),digits = 6)) %>% set_colnames(c("A1", "B1", "C1","A2", "B2", "C2","A3", "B3", "C3")),  
                  file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/full_zscores/seed", seedi, ".replicate",rep,".", i, ".gwas.csv"),
                  quote = FALSE,sep = ",", row.names = FALSE )
    }else
    {
      write.table(data.frame(round(cbind(first$se, second$se, third$se),digits = 6)) %>% set_colnames(c("A1", "B1", "C1","A2", "B2", "C2","A3", "B3", "C3")),  
                  file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/simulations_low_herit_sumstats/seed", seedi, ".replicate",rep,".", i, ".SE.csv"),
                  quote = FALSE,sep = ",", row.names = FALSE )
      write.table(data.frame(round(cbind(first$beta, second$beta, third$beta),digits = 6)) %>% set_colnames(c("A1", "B1", "C1","A2", "B2", "C2","A3", "B3", "C3")),  
                  file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/simulations_low_herit_sumstats/seed", seedi, ".replicate",rep,".", i, ".BETA.csv"),
                  quote = FALSE,sep = ",", row.names = FALSE )
    }
    
   
  }

}

