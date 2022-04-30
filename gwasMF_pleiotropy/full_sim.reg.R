source("/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/helper_functions_plieo.R")
pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr, cowplot, magrittr)
source("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
thou.g.mafs <- fread("/work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/plink.frq")
thou.g.pruned <- fread("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/500kb.0.04r2.prune.in", header = F)
thou.g.mafs.sub <- thou.g.mafs %>% filter(SNP %in% thou.g.pruned$V1)
thou.g.mafs.use <- thou.g.mafs.sub
rm(thou.g.mafs)


#Helper function
safeSVD <- function(pop1, pop2, pop3, nsnps)
{
  mat <- matrix(as.numeric(cbind(pop1, pop2, pop3)), nrow = nsnps)
  if(any(is.na(mat)))
  {
    message("NA in mat")
  }
  s <- as.matrix(data.frame(scale(mat)) %>% drop_na())
  if(dim(s)[1] != dim(mat)[1])
  {
    print("had to drop things....")
    print(i)
  }
  
  svd(s)
  
}









set.seed(10)

#FIRST: create the genotypes
message("creating genotypes for 3 populations...")
N_tot <- 50000
genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop1.geno <- t(do.call("cbind", genotypes.all))
rm(genotypes.all)

genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop2.geno <- t(do.call("cbind", genotypes.all))
rm(genotypes.all)


genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop3.geno <- t(do.call("cbind", genotypes.all))
rm(genotypes.all)


#Perform the test several times
#cycle through the active SNPs
n.snps <- 1000
for(i in 1:9)
{
  start.snp <- (n.snps*(i-1)) + 1
  end.snp <- n.snps * i
  print(start.snp)
  print(end.snp)

  overlap.list <- c(1,1000,5000,7500,10000,12500, 15000)
  pop1.null <- list()
  
  
  y1.null <- phenotypeBuilder(NULL, pop1.geno, 1, type = "0,0,0")
  pop2.null <- list() 
  y2.null <- phenotypeBuilder(NULL, pop2.geno, 1, type = "0,0,0")
  
  pop3.null <- list() #
  y3.null <- phenotypeBuilder(NULL, pop3.geno, 1, type = "0,0,0")
  n.fixed = 15000
  pca.each <- list()
  #start with 3 groups for simplicity in overlapping
  message("Creating gwas summary stats...")
  for(j in overlap.list)
  {
    print(j)
    pop1.null[[j]] <- relatedOverlapBETTER(n.snps, j, 3, y1.null, N.tot = 50000, 
                                           genotypes=pop1.geno[,start.snp:end.snp], fixed.count = n.fixed, parallel = 7)
    pop2.null[[j]] <- relatedOverlapBETTER(n.snps, j, 3, y2.null, N.tot = 50000, 
                                           genotypes=pop2.geno[,start.snp:end.snp], fixed.count = n.fixed,parallel = 7)
    
    pop3.null[[j]] <- relatedOverlapBETTER(n.snps, j, 3, y3.null, N.tot = 50000, 
                                           genotypes=pop3.geno[,start.snp:end.snp], fixed.count = n.fixed,parallel = 7)
    
    mat <- matrix(as.numeric(cbind(pop1.null[[j]], pop2.null[[j]], pop3.null[[j]])), nrow = n.snps)
    print(head(mat))
    if(any(is.na(mat)))
    {
      message("NA in mat")
    }
    s <- as.matrix(data.frame(scale(mat)) %>% drop_na())
    if(dim(s)[1] != dim(mat)[1])
    {
      print("had to drop things....")
      print(i)
    }
    
    #pca.each[[j]] <- svd(s)
    pca.each[[j]] <- safeSVD(pop1.null[[j]], pop2.null[[j]], pop3.null[[j]])
  }
  save(pca.each,file = paste0("/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/nullFAsim/3group_pca.", i, ".RData"))
  
  
  #Do the corresponding baseline runs with overlap of 0
  rand_indices.null <- c()
  rand_indices.adj.null <- c()
  pcas <- list()
  n.fixed <- 15000
  groups <- list()
  for(rep in 1:1000)
  {
    y1.baseline <- phenotypeBuilder(NULL, pop1.geno, 1, type = "0,0,0")
    y2.baseline <- phenotypeBuilder(NULL, pop2.geno, 1, type = "0,0,0")
    y3.baseline <- phenotypeBuilder(NULL, pop3.geno, 1, type = "0,0,0")
    n.snps <- 1000
    n.fixed = 15000
    #start with 2 groups for simplicity in overlapping
    pop1.n <- relatedOverlapBETTER(n.snps, 0, 3, y1.null, N.tot = 50000, 
                                   genotypes=pop1.geno[,start.snp:end.snp], fixed.count = n.fixed,parallel = 7)
    pop2.n <- relatedOverlapBETTER(n.snps, 0, 3, y2.null, N.tot = 50000, 
                                   genotypes=pop2.geno[,start.snp:end.snp], fixed.count = n.fixed, parallel = 7)
    pop2.n <- relatedOverlapBETTER(n.snps, 0, 3, y2.null, N.tot = 50000, 
                                   genotypes=pop3.geno[,start.snp:end.snp], fixed.count = n.fixed, parallel = )
    #evaluate
    mat <- matrix(as.numeric(cbind(pop1.n, pop2.n)), nrow = n.snps)
    pcas[[rep]] <- svd(scale(mat))
    groups[[rep]] <- cutree(tree = hclust(dist( pcas[[rep]]$v[,1])), k = 2) #use just one to do it.
    rand_indices.null <- c(rand_indices.null, rand.index(groups[[rep]], c(0,0,0,1,1,1)))
    rand_indices.adj.null <- c(rand_indices.adj.null, adj.rand.index(groups[[rep]], c(0,0,0,1,1,1)))
    
  }
}

if(FALSE)
{
  #Now, we read them in manually and do the things.
  #We also need to run the null one (no overlap)
  
    #now analyze all the results....
  library(fossil)
  message("FA on results...")
  k <- c(1,1,1,2,2,2,3,3,3)
  groups <- list()
  rand_index_real <- c()
  rand_index_real.adj <- c()
  for(j in seq_along(pca.each))
  {
    p <- pca.each[[j]]
    groups[[j]] <- cutree(tree = hclust(dist(p$v[,1])), k = 2)
    rand_index_real <- c(rand_index_real, rand.index(groups[[j]], k))
    rand_index_real.adj <- c(rand_index_real.adj, adj.rand.index(groups[[j]], k))
  }

#save these results out.....
  save(pca.each,rand_index_real,rand_index_real.adj, file = paste0("/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/nullFAsim/3group_pca.", i, ".RData"))
}