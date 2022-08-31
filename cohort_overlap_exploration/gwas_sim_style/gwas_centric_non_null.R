

set.seed(2647)
message("Seed set to 2647")
pacman::p_load(stats, data.table, magrittr, dplyr, tidyr, ggplot2, stringr, fossil)
source("/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/helper_functions_plieo.R")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
#load("/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/.RData")

#Generate genotypes for 2 populations...
thou.g.mafs <- fread("/scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/plink.frq")
thou.g.pruned <- fread("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/500kb.0.04r2.prune.in", header = F)
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



n.iter <- 50
rep.dat <- list()
message("Starting sim iters...")
n.snps <- 1000
snp.weights <- list()
snp.weights[[1]] <- c(rep(0.04,200), rep(0,800))
snp.weights[[2]] <- c(rep(0,200), rep(0.2,25),rep(0,775))
snp.weights[[3]] <- c(rep(0.05,10), rep(0,900),rep(0.07,90))
stopifnot(length(snp.weights[[3]]) == n.snps)
stopifnot(length(snp.weights[[2]]) == n.snps)
stopifnot(length(snp.weights[[1]]) == n.snps)
for(rep in 1:n.iter)
{
  start = Sys.time()  
  #Simulate new phenotypes
  overlap.list <- c(0,2250,4500,6750,9000,11250,13500,15000)
  pop1.null <- list()#genotypes.tab
  y1 <- phenotypeBuilder(snp.weights, pop1.geno[,1:n.snps], 1, type = "1+2=3")
  y2 <- phenotypeBuilder(snp.weights, pop2.geno[,1:n.snps], 1, type = "1+2=3")
  y3 <- phenotypeBuilder(snp.weights, pop3.geno[,1:n.snps], 1, type = "1+2=3")
  pop2.null <- list() #genotypes.tab.2
  pop3.null <- list()

  n.fixed = 15000
  #start with 2 groups for simplicity in overlapping
  for(i in overlap.list)
  {
    print(i)
    
    pop1.null[[i]] <- relatedOverlapBETTER(n.snps, i, 3, y1, N.tot = N_tot, 
                                           genotypes=pop1.geno[,1:n.snps], fixed.count = n.fixed)
    pop2.null[[i]] <- relatedOverlapBETTER(n.snps, i, 3, y2, N.tot = N_tot, 
                                           genotypes=pop2.geno[,1:n.snps], fixed.count = n.fixed)
    
    pop3.null[[i]] <- relatedOverlapBETTER(n.snps, i, 3, y3, N.tot = N_tot, 
                                           genotypes=pop3.geno[,1:n.snps], fixed.count = n.fixed)
  }
  
  pca.each <- lapply(overlap.list, function(i) svd(scale(cbind(pop1.null[[i]], pop2.null[[i]],pop3.null[[i]]))))
  groups <- list()
  rand_index_real <- c()
  rand_index_real.adj <- c()
  cor.def <- NULL
  for(i in seq_along(pca.each))
  {
    p <- pca.each[[i]]
    pve <- p$d^2/sum(p$d^2)
    write.table(data.frame(pve), 
                file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/",rep,".", overlap.list[i], ".pve.csv"),quote = FALSE,sep = ",", row.names = FALSE )
    write.table(data.frame(p$u), 
                file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/",rep,".", overlap.list[i], ".u.csv"),quote = FALSE,sep = ",", row.names = FALSE )
    write.table(data.frame(p$v), 
                file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/",rep,".", overlap.list[i], ".v.csv"),quote = FALSE,sep = ",", row.names = FALSE )
  }

  message("Itertion took ",  round(difftime(Sys.time(), start, units = "mins"), digits = 3), " mins")
}

