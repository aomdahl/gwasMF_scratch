################################################################################################################################
## Ashton Omdahl, November 2024

##
################################################################################################################################
pacman::p_load(magrittr, dplyr, ggplot2, data.table, ggsignif)


library(devtools)
library(optparse)
 `%>%` <- magrittr::`%>%`
writeOutResults <- function(V,U,traits, snps,outdir)
{
  if(length(ncol(V)) == 0) #If we get nothing back
  {
    file.create( paste0(outdir, "latent.factors.txt"))
    file.create(paste0(outdir, "latent.loadings.txt"))
  }else
  {
    V.df <- data.frame(traits, V) %>% magrittr::set_colnames(c("Study", paste0("V", 1:ncol(V))))
    write.table(V.df, file = paste0(outdir, "latent.factors.txt"), quote = FALSE, row.names = FALSE, sep = " ")

    U.df <- data.frame(snps, U) %>% magrittr::set_colnames(c("SNP", paste0("U", 1:ncol(V))))
    write.table(U.df, file = paste0(outdir, "latent.loadings.txt"), quote = FALSE, row.names = FALSE, sep = " ")

    write.table(traits, file =paste0(outdir, "trait_out_order.txt"), quote = FALSE, row.names = FALSE,col.names = FALSE)
  }

}



option_list <- list(
make_option(c("--gwas_effects"), type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait"),
make_option(c("--uncertainty"), type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait"),
make_option(c("--lambda_gc"), type = 'character', help = "Specify the path to the genomic correction coefficients. If none provided, none used", default = ""),
make_option(c("--trait_names"), type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the order in the input tables.", default = ""),
make_option(c("--weighting_scheme"), type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE"),
make_option(c("--covar_matrix"), type = 'character', help = "Path to LDSC estimates of covariance effect. No adjustment made if none provided", default = ""),
make_option(c("-c", "--converged_obj_change"), type = 'numeric', help = "Specify the objective percent change required to achieve convergence", default = 0.001),
make_option(c("-i", "--niter"), type = 'numeric', help = "Cap the number of iterations on the matrix learning step", default = 300),
make_option(c("--outdir"), type = "character", help = "Source file location"),
make_option(c("--drop_phenos"), type = "character", help = "Specify phenotypes to exclude (useful if encountering issues with covariance adjustment)", default = ""),
make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
make_option(c("--debug"), type = "logical", help = "if want debug run", action = "store_true", default = FALSE),
make_option(c("--overview_plots"), type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE),
make_option(c("-K", "--nfactors"), type = "character", help = "specify the number of factors. Options are a number, or MAX, KAISER, K-2, CG", default = "GRID"),
make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
make_option(c("--bic_var"), type = 'character', help = "Specify the bic method to use. Options are [sklearn_eBIC,sklearn,Zou_eBIC,Zou,dev_eBIC,dev, std,NONE]", default = "sklearn_eBIC"),
make_option(c("-p", "--param_conv_criteria"), type = 'character', help = "Specify the convergene criteria for parameter selection", default = "BIC.change"),
make_option(c("-r", "--rg_ref"), type = 'character', help = "Specify a matrix of estimated genetic correlations to initialize V from", default = ""),
make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud"),
make_option(c("-n", "--ncores"), type="integer", default= 1, help="Do you want to run on multiple cores? Only influences the K initialization step at this point."),
#These related to covariance matrix tweaks
make_option(c("-s", "--sample_sd"), type="character", default= "", help="File containing the standard deviation of SNPs; if provided, used to scale LDSC gcov terms."),
make_option(c("-b", "--block_covar"), type = "numeric", default= 0.2, help="Specify the degree to which a block structure is enforced, by cluster distance. Default is 0.2"),
make_option(c("--covar_se_matrix"), type = "character", default="", help="Path to covar se matrix, needed if using Strimmer gamma specification."),
make_option(c("--intermediate_bic"), type = "character", default="", help="OPTIONAL:Path to intermediate BIC file, will initiate factorization from here. "),
make_option(c("--subset_seed"), type = "numeric", default=-1, help="Specify a seed to subset the data for a variant test. For debugging purposes."),
make_option(c("--model_select"), type = "logical", default=FALSE,action ="store_true", help="Specify this if you only wish to do the model selection step."),
make_option(c("-g", "--WLgamma"), type="character", default= "0",
            help="Specify the extent of the WL shrinkage on the input covariance matrix.\nCan pass as a number (1 for no covariance adjustment, 0 for full adjustment), or to specify a method (either MLE or Strimmer)")
)



t=c("--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.se.tsv",
    "--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.beta.tsv",
"--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete.trait_list.tsv",
    "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//ldsr_results/panUKBB_complete/summary_data/gcov_int.tab.csv",
"--outdir=/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/",
"--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/sample_sd_report.tsv","--fixed_first",
"--WLgamma=Strimmer","--covar_se_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv",
"--nfactors=134", "--bic_var=sklearn_eBIC", "--converged_obj_change=0.001", "--WLgamma=Strimmer", "--model_select")


#args_input <- parse_args(OptionParser(option_list=option_list))#, args = t)
args_input <- parse_args(OptionParser(option_list=option_list), args = t)
#TODO: cut oout this middle-man, make a more efficient argument input vector.

#setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/"); load_all()
library(gleanr)
args <-fillDefaultSettings(args_input)

writeRunReport(args)

#Establish the internal settings
option <- readInSettings(args)
option$alpha1 <- NA
option$lambda1 <- NA
outdir <-args$outdir

input.dat <- readInData(args)
X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c; C <- input.dat$C_block
option$C <- C
option$WLgamma <- args$WLgamma


#We need both the final V and the sparsity selected
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
#Prep regression datasets (weights, x, etc.)
#Keep these things in memory so not recreating each time
reg.vect <- prepRegressionElements(X,W,W_c,option) #returns sparse matrices since that saves a ton of memory by doing this once up front
#Loading that all in each time was slow...
userMessage(args$verbosity, "Outcome data weighted and transformed for regression", thresh=0)
opath=outdir

## Shuffle V
set.seed(2)
shuffV <- ret$V[sample(1:nrow(ret$V), nrow(ret$V)),]
## Not sure which way, for now let's just be conservative and shuffle all jointly

#What I should do is initialize a GLEANER object, and have all these things in there....
#Can I make the regvect elements sparse?
#Perform sparsity selection,

#Transition over to the full run:
message("Convergence set to ",  option$conv0 )
option$K <- ncol(ret$V)
option$alpha1 <- ret$autofit_alpha[1]#This is applied to U
option$lambda1 <- ret$autofit_lambda[1] #This is applied to V
option$iter <- 0.5 #only fit U

#Perform optimization
ret.permute <- gwasML_ALS_Routine(option, X, W, W_c, shuffV, reg.elements=reg.vect, no.prune = TRUE) #I like this better
image(ret.permute$U)
plotCorrelationHeatmap(cor(ret.permute$U, ret$U))
plotCorrelationHeatmap(cor(ret.permute$U, ret$U))
plotCorrelationDiagMatch(ret$U,ret.permute$U)
plotCorrelationDiagMatch(ret$V,shuffV) + ggtitle("Original V vs shuffled V")

#A more aggressive shuffling:
set.seed(3)
shuffV.alt <- apply(ret$V, 2, function(x)sample(x,size=length(x), replace=FALSE))
plotCorrelationDiagMatch(ret$V,shuffV.alt) + ggtitle("Original V vs column-shuffled V")
ret.permute.agg <- gwasML_ALS_Routine(option, X, W, W_c, shuffV.alt, reg.elements=reg.vect, no.prune = TRUE) #I like this better
plotCorrelationDiagMatch(ret$U,ret.permute.agg$U) + ggtitle("Original U vs alt shuffled U")

#So I would just need to repeat this 100x and see what happens?
#hmmm very slow.
shuff.u <- list()
for(i in 1:50)
{
  print(i)
  set.seed(i)
  shuffV.alt <- apply(ret$V, 2, function(x)sample(x,size=length(x), replace=FALSE))
  shuff.u[[i]] <- gwasML_ALS_Routine(option, X, W, W_c, shuffV.alt, reg.elements=reg.vect, no.prune = TRUE) #I like this better
}
save(shuff.u, file="/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/permute_50.RData")
#load("/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/permute_50.RData")
#import functinos from getTopGenes
## Now get the p-values:
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")
source("/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/src/get_pLI.R")

#Get blood traits out
blood.trait.idx <- which(ret$V[,4] != 0)
blood.traits <- ret$trait.names[blood.trait.idx]
blood.idx <- which(ret$trait.names %in% blood.traits)
blood.factors <- unique(c(apply(ret$V[blood.idx,], 1, function(x) which(x != 0))))

#I did it once and lost it. dummy.
## Now for the enrichment test:
#JUST THE IDP GENES
early_Megakaryopoiesis=c("THPO","MPL","GATA1","RUNX1","FLI1","ETV6","GFI1B","HOXA11","MECOM","ANKRD26","RBM8A")
late_mk <- c("AP3B1","HPS1","HPS3","HPS4","HPS5","HPS6","BLOC1S3","BLOC1S6","DTNBP1","LYST","VPS33B","VIPAS39","STXBP2", "NBEA", "NBEAL2", "CYCS","SRC","SLFN14","PLAU","STIM1")
protoplatelet <- c("MYH9", "WAS", "ACTN1", "FLNA", "TUBB1", "DIAPH1", "GP1BA", "GP1BB","GP9","ITGA2B", "ITGB3", "VWF")
platelet_function <- c("P2RY12","TBXA2R", "TBXAS1", "PLA2G4A"," ITGA2B"," RASGRP2","VWF", "ITGB3","FERMTS", "GP1BA","GP9","GP1BB","GP6"," ANO6")
macro.genes <- c("THPO","ANKRD26","ETV6","FLI1","GATA1","GFI1B","HOXA11","MECOM","RUNX1","NBEAL2","ABCG5","ABCG8","GNE","MPIG6B","SLFN14","SRC","ACTB","ACTN1","CDC42","DIAPH1","FLNA","MYH9","TUBB1","GP1BA","GP1BB","GP9","ITGA2B","ITGB3","VWF")

#ADDING IN THE ONES RELATED TO HEREDITARY THROMBOCYTOPENIA (WHERE MISSING)
#OMITTING the grey ones with rare variants causing hereditary HT
em_hmt <- c()
mk_hmt <- c("ABCG5","ABCG8","GNE","MPIG6B")
pp_hmt <- c("ACTB","CDC42")
plt_hmt <- c()


f4.tp <- loadGeneEnrichments(4, "top_fe")
f11.tp <- loadGeneEnrichments(11, "top_fe")
f23.tp <- loadGeneEnrichments(23, "top_fe")
full.bg <- unique(f4.tp$bg_genes)
#For viz- the scree plots of each?
#prioritize_snps from get top Genes
#Need to prioritize SNPs
#only care about 3 factors here, so let's look at those
f <- c(4,11,23)
factor=11
test_sets <- list(macro.genes,early_Megakaryopoiesis,c(late_mk,mk_hmt), c(protoplatelet,pp_hmt),platelet_function)
null.f4 <- lapply(test_sets, function(x) list())
null.f11 <- lapply(test_sets, function(x) list())
for(i in 1:50)
{
  top.snps <- prioritize_snps(as.matrix(shuff.u[[i]][,factor]),snp_ids = ret$snp.ids,method="top_fe")
  for(j in 1:length(test_sets))
  {
    #null.f4[[j]][[i]]  <- testEnrichment(unique(top.snps[[1]]), unique(full.bg), test_sets[[j]],conf.level=0.90,alternative="two.sided")
    null.f11[[j]][[i]]  <- testEnrichment(unique(top.snps[[1]]), unique(full.bg), test_sets[[j]],conf.level=0.90,alternative="two.sided")
  }

}
#The moment of truth- distribution of hypergeometric p-values
p.vals. <- sapply(null.f4, function(x) x$test$p.value)
ors.all <- sapply(null.f4, function(x) x$test$estimate)
hist(sapply(null.f4[[1]], function(x) x$test$p.value),main="Null P, HMTP genes")
hist(sapply(null.f4[[2]], function(x) x$test$p.value),main="Null P, Early Mega. genes")
hist(sapply(null.f4[[3]], function(x) x$test$p.value),main="Null P, Late Mega. genes")
hist(sapply(null.f4[[4]], function(x) x$test$p.value),main="Null P, Protoplatelet formation genes")
hist(sapply(null.f4[[5]], function(x) x$test$p.value),main="Null P, Platelt genes")

hist(sapply(null.f11[[1]], function(x) x$test$p.value),main="Null P, HMTP genes")
hist(sapply(null.f11[[2]], function(x) x$test$p.value),main="Null P, Early Mega. genes")
hist(sapply(null.f11[[3]], function(x) x$test$p.value),main="Null P, Late Mega. genes")
hist(sapply(null.f11[[4]], function(x) x$test$p.value),main="Null P, Protoplatelet formation genes")
hist(sapply(null.f11[[5]], function(x) x$test$p.value),main="Null P, Platelt genes")



#Really really really no evidence of inflation here. Nothing is coming up at all.
hist(ors.all)
#okay, really really nothing going on here. Cool.
f4.tp$all_genes %>% dplyr::select(RSID, gene) %>% mutate("top_snp"=RSID %in% top.snps[[1]]) %>% distinct() %>% arrange(top_snp)

#Permuted and preferred test.

#Better names:
factor_sets = list(f4.tp$set_genes, f11.tp$set_genes, f23.tp$set_genes)
factor_titles <- c("F4","F11","F23")
or_test_hmt <- testEnrichment(unique(f4.tp$set_genes), unique(full.bg), macro.genes,conf.level=0.90,alternative="two.sided")
plot.me <- data.frame("OR"=or_test_hmt$test$estimate, "upper"=or_test_hmt$test$conf.int[2], "lower"=or_test_hmt$test$conf.int[1], "p"=or_test_hmt$test$p.value, "Factor"="U4")
#Quick nice plot of this 10/28
ggplot(plot.me, aes(x=Factor, y=OR)) + geom_point(size=3) +
  geom_errorbar( aes(x=Factor, ymin=lower, ymax=upper), width=0.3, colour="black", alpha=0.9) + theme_classic(16) +
  geom_hline(yintercept = 1,color="black", lty="dashed")+ ylab("Enrichment OR")  + coord_flip()

#Sandbox- could show as a z-score
#qnorm(or_test_hmt$test$p.value*2,lower.tail = TRUE)

##Harder test- just within blood genes
blood.bg <- c(f4.tp$set_genes,f11.tp$set_genes,f23.tp$set_genes)
testEnrichment(unique(f4.tp$set_genes), unique(blood.bg), macro.genes,conf.level=0.90,alternative="two.sided")
testEnrichment(unique(f11.tp$set_genes), unique(blood.bg), macro.genes,conf.level=0.90,alternative="two.sided")
testEnrichment(unique(f23.tp$set_genes), unique(blood.bg), macro.genes,conf.level=0.90,alternative="two.sided")



## Next test- all IDP genes at stages
phases=c("Early megakaryopoiesis","Late megakaryopoiesis", "Protoplatelet formation", "Platelet function")
test_sets <- list(early_Megakaryopoiesis,c(late_mk,mk_hmt), c(protoplatelet,pp_hmt),platelet_function)
#tabEnrichmentTestPermuted <- function(factor_sets,full.bg,factor_titles, sets_to_test,gene_set_names,n_permute=1000)
if(LOAD)
{
  load("/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/enrichment_test_platelets.RData")
  
}else
{
  ##First test- just those for macro genes
  permuted.joint.hmtc <- tabEnrichmentTestPermuted(factor_sets,full.bg,factor_titles, list(macro.genes),c("HMTC"),n_permute=10000,conf.level=0.9)
  permuted.joint.bg <- tabEnrichmentTestPermuted(factor_sets,full.bg,factor_titles, test_sets,phases,n_permute=50000,conf.level=0.9)
  save(permuted.joint.bg,permuted.joint.hmtc,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/enrichment_test_platelets.RData")
  
}


permuted.joint.bg$Phase <- factor(permuted.joint.bg$Phase, levels=phases)
permuted.joint.bg$Factor <- factor(permuted.joint.bg$Factor, levels = paste0("F",1:100))
permuted.joint.bg$fdr_sig <- ifelse(permuted.joint.bg$FDR < 0.05, "FDR < 0.05", "FDR > 0.05")



