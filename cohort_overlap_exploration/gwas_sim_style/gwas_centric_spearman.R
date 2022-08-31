set.seed(2648)
message("Seed set to 2647 normally, moved for this run.")
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
#genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
#pop3.geno <- t(do.call("cbind", genotypes.all))
rm(genotypes.all)



n.iter <- 42
rep.dat <- list()
message("Starting sim iters...")
for(rep in 1:n.iter)
{
  start = Sys.time()  
  #Simulate new phenotypes
    overlap.list <- c(1,1000,5000,7500,10000,12500)
    pop1.null <- list()#genotypes.tab
    y1.null.similar <- phenotypeBuilder(NULL, pop1.geno, 1, type = "1+2=3")
    y2.null.similar <- phenotypeBuilder(NULL, pop2.geno, 1, type = "1+2=3")
    #y3.null.similar <- phenotypeBuilder(NULL, pop2.geno, 1, type = "1+2=3")
    pop2.null <- list() #genotypes.tab.2
    #pop3.null <- list()
    n.snps <- 1000
    n.fixed = 15000
    #start with 2 groups for simplicity in overlapping
    for(i in overlap.list)
    {
      print(i)

      pop1.null[[i]] <- relatedOverlapBETTER(n.snps, i, 3, y1.null.similar, N.tot = N_tot, 
                                             genotypes=pop1.geno, fixed.count = n.fixed)
      pop2.null[[i]] <- relatedOverlapBETTER(n.snps, i, 3, y2.null.similar, N.tot = N_tot, 
                                             genotypes=pop2.geno, fixed.count = n.fixed)
      
      #pop3.null[[i]] <- relatedOverlapBETTER(n.snps, i, 3, y2.null.similar, N.tot = N_tot, 
      #                                       genotypes=pop3.geno, fixed.count = n.fixed)
    }

    pca.each <- lapply(overlap.list, function(i) svd(scale(cbind(pop1.null[[i]], pop2.null[[i]]))))
    groups <- list()
    rand_index_real <- c()
    rand_index_real.adj <- c()
    cor.def <- NULL
    for(i in seq_along(pca.each))
    {
      p <- pca.each[[i]]
      pve <- p$d^2/sum(p$d^2)
      rownames(p$v) <- c("A1", "B1", "C1", "A2", "B2", "C2")
      y=c(0, 0, 0, 1, 1, 1)
      sig.pves <- which(pve > (1/length(p$d)))
      corlook.cor <- lapply(sig.pves, function(x) cor.test(y=y, x = p$v[,x], method = "spearman")$estimate)
      corlook.p <- lapply(sig.pves, function(x) cor.test(y=y, x = p$v[,x])$p.value)
      cor.def <- rbind(cor.def, data.frame("rho" = unlist(corlook.cor), "p"= unlist(corlook.p), "pc" = sig.pves, "overlap" = overlap.list[i]))
        write.table(data.frame(pve),
                file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/",rep,".", overlap.list[i], ".pve.csv"),quote = FALSE,sep = ",", row.names = FALSE )
    write.table(data.frame(p$u),
                file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/",rep,".", overlap.list[i], ".u.csv"),quote = FALSE,sep = ",", row.names = FALSE )
    write.table(data.frame(p$v),
                file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/",rep,".", overlap.list[i], ".v.csv"),quote = FALSE,sep = ",", row.names = FALSE )
  }

    }
    write.table(cor.def, file = paste0("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/iter",rep,"res.csv"),quote = FALSE,sep = ",", row.names = FALSE )
    
    #ggplot(cor.def, aes(y = abs(rho), x = as.factor(round((overlap/150), digits =3)), fill =  as.factor(pc))) + 
    #  geom_bar(stat = "identity", position = "dodge") + theme_classic(15)  + 
    #  labs(fill = "Which PC") + ggtitle("Top source of variation dominated by cohort overlap as\nsample overlap increases") + xlab("Percent overlap in cohort")
  message("iter ",rep," complete")
  message("Itertion took ",  round(difftime(Sys.time(), start, units = "mins"), digits = 3), " mins")
}

