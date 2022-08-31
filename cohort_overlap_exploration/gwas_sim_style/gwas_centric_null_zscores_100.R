args = commandArgs(trailingOnly=TRUE)
set.seed(args[1])
seedi = args[1]
message("Seed set to ", args)
message("this version of the script simply saves out all the GWAs summary stats, doesn't actually do any analysis.")
pacman::p_load(stats, data.table, magrittr, dplyr, tidyr, ggplot2, stringr, fossil)
source("/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/helper_functions_plieo.R")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
#load("/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/.RData")

#Generate genotypes for 3 populations...
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

write.table(thou.g.mafs.use, file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/full_zscores/seed", seedi, "MAF.csv"), quote =FALSE,sep = ",", row.names = FALSE )

n.iter <- 100
rep.dat <- list()
message("Starting sim iters...")
n.snps <- 1000
for(rep in 1:n.iter)
{
  start = Sys.time()  
  #Simulate new phenotypes
  overlap.list <- c(0,1000,5000,7500,10000,12500)
  pop1.null <- list()#genotypes.tab
  y1 <- phenotypeBuilder(NULL, pop1.geno, 1, type = "1+2=3")
  y2 <- phenotypeBuilder(NULL, pop2.geno, 1, type = "1+2=3")
  pop2.null <- list() #genotypes.tab.2

  n.fixed = 15000
  #start with 2 groups for simplicity in overlapping
  for(i in overlap.list)
  {
    print(i)
    
    first <- relatedOverlapBETTER(n.snps, i, 3, y1, N.tot = N_tot, 
                                           genotypes=pop1.geno[,1:n.snps], fixed.count = n.fixed)

   second <- relatedOverlapBETTER(n.snps, i, 3, y2, N.tot = N_tot, 
                                           genotypes=pop2.geno[,1:n.snps], fixed.count = n.fixed)
    

    write.table(data.frame(round(cbind(first, second),digits = 6)) %>% set_colnames(c("A1", "B1", "C1","A2", "B2", "C2")),  
                file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/full_zscores/seed", seedi, ".replicate",rep,".", i, ".gwas.csv"),
                quote = FALSE,sep = ",", row.names = FALSE )
  }
  message("Finished iteration ", n.iter)

}

