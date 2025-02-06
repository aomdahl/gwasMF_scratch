LOAD=TRUE
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")
source("/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/src/get_pLI.R")

pacman::p_load(magrittr, dplyr, ggplot2, data.table, ggsignif)


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



load("/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed23_/permute_100.RData")
#get SNP order:
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
v2g <- load_v2G() %>% mutate("rsid" = factor(rsid, levels = ret$snp.ids)) %>% arrange(rsid) #then force the order, then get the gene id then do your test
stopifnot(!any(is.na(v2g$rsid)))
shuff_results_U[[1]]
#This seems like it works.
#We have our procedure for dsignating 
testEnrichmentLogit(platelet_function,v2g$gene_id, abs(ret$U[,4]), collapse = "max")
testEnrichment(unique(f4.tp$set_genes), unique(f4.tp$bg_genes), platelet_function,conf.level=0.90,alternative="two.sided")

#only care about 3 factors here, so let's look at those
f <- c(4,11,23)
factor=11
test_sets <- list(macro.genes,early_Megakaryopoiesis,c(late_mk,mk_hmt), c(protoplatelet,pp_hmt),platelet_function)
null.f4 <- lapply(test_sets, function(x) list())
null.f11 <- lapply(test_sets, function(x) list())
null.f11.logit.max <- lapply(test_sets, function(x) list())
null.f11.logit.mean <- lapply(test_sets, function(x) list())
for(i in 1:length(shuff_results_U))
{
  top.snps <- prioritize_snps(as.matrix(shuff_results_U[[i]][,factor]),snp_ids = ret$snp.ids,method="top_fe")
  for(j in 1:length(test_sets))
  {
    #null.f4[[j]][[i]]  <- testEnrichment(unique(top.snps[[1]]), unique(full.bg), test_sets[[j]],conf.level=0.90,alternative="two.sided")
    null.f11[[j]][[i]]  <- testEnrichment(unique(top.snps[[1]]), unique(full.bg), test_sets[[j]],conf.level=0.90,alternative="two.sided")
    null.f11.logit.max[[j]][[i]] <- testEnrichmentLogit(test_sets[[j]],v2g$gene_id, abs(shuff_results_U[[i]][,factor]), collapse = "max")
    null.f11.logit.mean[[j]][[i]] <- testEnrichmentLogit(test_sets[[j]],v2g$gene_id, abs(shuff_results_U[[i]][,factor]), collapse = "mean")
  }
  
}
#Now parse a given matrix or set thereof
j=1
permute.dat <- cbind(
do.call("rbind", lapply(null.f11[[j]], function(x)
  {
  c(x$test$p.value,x$test$estimate,x$logOR,x$logOR_SE)
}))[1:46,],
do.call("rbind", lapply(null.f11.logit.max[[j]], function(x)
{
  c(summary(x$test)$coef[8],x$logOR,x$logOR_SE)
})),
do.call("rbind", lapply(null.f11.logit.mean[[j]], function(x)
{
  c(summary(x$test)$coef[8],x$logOR,x$logOR_SE)
}))
)
colnames(permute.dat) <- c("fisher_p", "fisher_OR", "fisher_logOR", "fisher_logOR_SE", 
                           "maxlogit_p", "maxlogit_b", "maxlogit_SE",
                           "avglogit_p", "avglogit_b", "avglogit_SE")
source("/data/abattle4/aomdahl1/gene_utils/qqplot.R")
qqunif.plot(permute.dat$maxlogit_p)
qqunif.plot(permute.dat$avglogit_p)
qqunif.plot(permute.dat$fisher_p)
#erm that's not encouraging
#logit may be inflated, fisher is overly conservative.


null.f11[[j]][[1]]$test$p.value
null.f11[[j]][[1]]$test$estimate
null.f11[[j]][[1]]$logOR
null.f11[[j]][[1]]$logOR_SE
testEnrichmentLogit(platelet_function,v2g$gene_id, abs(shuff_results_U[[1]][,4]), collapse = "max")
  
  
  
  
f4.tp <- loadGeneEnrichments(4, "top_fe")
f11.tp <- loadGeneEnrichments(11, "top_fe")
f23.tp <- loadGeneEnrichments(23, "top_fe")


#Keep the background consistent with published tests, don't want to look shdy.
nice.background <- unique(c((f4.tp$all_genes %>% dplyr::filter(U != 0))$RSID,(f23.tp$all_genes %>%  dplyr::filter(U!=0))$RSID), (f11.tp$all_genes %>%  dplyr::filter(U!=0))$RSID)
full.bg <- (f4.tp$all_genes %>%  dplyr::filter(RSID %in% nice.background) %>% group_by(RSID) %>% slice_head(n=1) %>% ungroup() %>%  dplyr::filter(!is.na(gene)))$gene

f4.all <- testEnrichment(unique(f4.tp$set_genes), unique(full.bg), macro.genes,conf.level=0.90,alternative="two.sided")
f11.all <- testEnrichment(unique(f11.tp$set_genes), unique(full.bg), macro.genes,conf.level=0.90,alternative="two.sided")
f23.all <- testEnrichment(unique(f23.tp$set_genes), unique(full.bg), macro.genes,conf.level=0.90,alternative="two.sided")


#A background alternative we could try
tops <- f4.tp$gene_snp_singular %>% select(RSID, gene, U, rank, top_snp) 
top.genes <- unique((tops %>% filter(top_snp))$gene)
all(top.genes %in% f4.tp$set_genes)
not.top.genes <- unique((tops %>% filter(!top_snp))$gene)
testEnrichment(unique(f4.tp$set_genes), unique(not.top.genes), macro.genes,conf.level=0.90,alternative="two.sided")
#This emphasizes the enrichment even more, because some genes get pulled from the BG.
calcOR(f4.all$tab)
f4.all
#I am being conservative. Might matter in some cases.
#LOOKING AT NEW ENRICHMENT STUFF.
factor_sets = list(f4.tp$set_genes, f11.tp$set_genes, f23.tp$set_genes)
factor_titles <- c("F4","F11","F23")
phases=c("Early megakaryopoiesis","Late megakaryopoiesis", "Protoplatelet formation", "Platelet function")
test_sets <- list(early_Megakaryopoiesis,c(late_mk,mk_hmt), c(protoplatelet,pp_hmt),platelet_function)

  all.tests <- list()
  customOR <- list()
  j=1
  for(i in 1:length(factor_sets))
  {
    message("Factor set ",i, " of ", length(factor_sets) )
    f <- factor_sets[[i]]
    fname <- factor_titles[[i]]
    all.tests[[factor_titles[[i]]]] <- lapply(test_sets, function(set)  testEnrichment(unique(f), unique(full.bg), set,conf.level=0.90,alternative="two.sided"))
    for(test in all.tests[[factor_titles[[i]]]])
    {
      customOR[[j]] <- calcOR(test$tab)
      j=j+1
    }
  }

  sapply(all.tests, function(x) sapply(x, function(y) y$test$estimate))
Or.all <- c(sapply(all.tests, function(x) sapply(x, function(y) y$test$estimate)))
  
custom.version <- sapply(customOR, function(x)x$OR) 
plot(Or.all, custom.version); abline(b=1,a=0, col="blue") 
#Perfect. Now.
custom.version.beta <- sapply(customOR, function(x)x$log_OR) 
custom.version.se <- sapply(customOR, function(x)x$SE) 
test.ash <- ashr::ash(custom.version.beta,custom.version.se)
test.ash$result
get_lfsr(test.ash)

#logistic regression enrichment version- testing
or_test_hmt <- testEnrichment(unique(f4.tp$set_genes), unique(full.bg), macro.genes,conf.level=0.90,alternative="two.sided")
in.set.genes <- as.numeric(unique(full.bg) %in% unique(f4.tp$set_genes))
test.set.genes <- as.numeric(unique(full.bg) %in%macro.genes)
## Approach 1) 
summary(glm(test.set.genes~in.set.genes,family = "binomial"))
max.effect.by.gene <- f4.tp$all_genes %>% group_by(gene) %>% 
  slice_max(abs(U)) %>% mutate("in_set" = as.numeric(gene %in% macro.genes)) %>% select(gene, U, in_set) %>% distinct()

summary(glm(in_set ~ abs(U),family="binomial",data=max.effect.by.gene))


#So we get the local false sign rate of the log-odds ratio?

calcOR <- function(data,alpha=0.05)
{
  OR=(data[1]*data[4])/(data[2]*data[3])
  SE=sqrt(sum(1/as.vector(data)))
  list( "OR"=OR, "SE"=SE, "log_OR"=log(OR),
  "logOR_CI"= log(OR) + c(qnorm(1-(alpha/2)) * SE,qnorm(1-(alpha/2)) * SE*-1))
}


######### What does it look like for Panglao?

###Panglao DB:
all.markers <- fread("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
#unique(all.markers$`cell type`)
platelet.genes <- (dplyr::filter(all.markers, `cell type` == "Platelets"))$`official gene symbol`
hsc.genes <- (dplyr::filter(all.markers, `cell type` == "Hematopoietic stem cells"))$`official gene symbol`
mk.genes <- (dplyr::filter(all.markers, `cell type` == "Megakaryocytes"))$`official gene symbol`

factor_sets = list(f4.tp$set_genes, f11.tp$set_genes, f23.tp$set_genes)
factor_titles <- c("F4","F11","F23")
sets_to_test <- list(platelet.genes,hsc.genes,mk.genes)
gene_set_names <- c("Platelets","HSCs", "Megakaryocytes")

unique.permuted.cell.types <- tabEnrichmentTestPermuted(factor_sets,f11.tp$bg_genes,factor_titles, sets_to_test,gene_set_names,n_permute=10000,conf.level=0.9)


all.tests.panglao <- list()
customOR.panglao <- list()
j=1
full.bg <- f11.tp$bg_genes
for(i in 1:length(factor_sets))
{
  message("Factor set ",i, " of ", length(factor_sets) )
  f <- factor_sets[[i]]
  fname <- factor_titles[[i]]
  all.tests.panglao[[factor_titles[[i]]]] <- lapply(sets_to_test, function(set)  testEnrichment(unique(f), unique(full.bg), set,conf.level=0.90,alternative="two.sided"))
  for(test in all.tests.panglao[[factor_titles[[i]]]])
  {
    customOR.panglao[[j]] <- calcOR(test$tab)
    j=j+1
  }
}

#Visualize:
sapply(all.tests.panglao, function(x) sapply(x, function(y) y$test$estimate))
Or.all.panglao <- c(sapply(all.tests.panglao, function(x) sapply(x, function(y) y$test$estimate)))

custom.version.panglao <- sapply(customOR.panglao, function(x)x$OR) 
plot(Or.all.panglao, custom.version.panglao); abline(b=1,a=0, col="blue")  #looks good

custom.version.beta.pang <- sapply(customOR.panglao, function(x)x$log_OR) 
custom.version.se.pang <- sapply(customOR.panglao, function(x)x$SE) 
test.ash.pang <- ashr::ash(custom.version.beta.pang,custom.version.se.pang)
test.ash.pang$result
get_lfsr(test.ash.pang)
#ASSIGNMENTS ARE WRONG_ revist
make_dat <- data.frame(expand.grid(c("HSC","Megakaryocyte", "Platelets"),c("F4","F11","F23"))) %>% mutate("lfsr"=get_lfsr(test.ash.pang), "OR"=Or.all.panglao)
barplot(get_lfsr(test.ash.pang))
ggplot(make_dat, aes(x=Var2, y=lfsr)) + geom_bar(stat="identity") + facet_wrap(~Var1)

ggplot(make_dat, aes(x=Var2, y=OR)) + geom_bar(stat="identity") + facet_wrap(~Var1)
#I think this is wrong, need t



#####
real.deal <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/assess_test1/enrichment_results.tsv")
ash.look.max.logit <- ash(real.deal$maxlogit_b, real.deal$maxlogit_SE)
perf.max.logit <- ash.look.max.logit$result %>% mutate("factor"=real.deal$factor, "gene_set"=real.deal$gene_set)

ash.look.avg.logit <- ash(real.deal$avglogit_b, real.deal$avglogit_SE)
perf.avg.logit <- ash.look.avg.logit$result %>% mutate("factor"=real.deal$factor, "gene_set"=real.deal$gene_set)

ash.look.fisher.or <- ash(log(real.deal$fisher_OR + 0.001), real.deal$fisher_logOR_SE)
perf.fisher.or <- sh.look.fisher.or$result %>% mutate("factor"=real.deal$factor, "gene_set"=real.deal$gene_set)




cowplot::plot_grid(
  ggplot(perf.max.logit,aes(x=as.factor(factor), y=PositiveProb, fill=lfsr)) + geom_bar(stat="identity") + facet_wrap(~gene_set, nrow=1) + theme_bw() + ggtitle("max logit"),
  ggplot(perf.avg.logit,aes(x=as.factor(factor), y=PositiveProb, fill=lfsr)) + geom_bar(stat="identity") + facet_wrap(~gene_set, nrow=1) + theme_bw() + ggtitle("avg logit"),
  ggplot(perf.fisher.or,aes(x=as.factor(factor), y=PositiveProb, fill=lfsr)) + geom_bar(stat="identity") + facet_wrap(~gene_set, nrow=1) + theme_bw() + ggtitle("fisher OR"),nrow=3)

### need to think about which background to use
### most generous: all SNPs included in analysis
## more restrictive- those that are zeroed out in that particular factor
## but that is part of the procedure we are testing, so maybe that's not fair.
## The complement of the genes in the set- so all those genes not mapped to be a top SNP
  ## This is different, but does account for teh fact that some genes could be in both sets.
  ## This strikes me as a more approriate way to go about this.
## Does it eben matter though?

  


###Looking at the null.
permuted <- fread("/scratch16/abattle4/ashton/snp_networks//scratch/manuscript_reviews/permute_testing/assess_test1/permuted_enrichment_results.tsv")
View(permuted)
##
ggplot(permuted, aes(x=fisher_p)) + geom_histogram() + 
  facet_grid(factor~gene_set) + theme_bw() + ggtitle("Fisher P-values")
#so no issue with inflation to be sure.....
ggplot(permuted, aes(x=avglogit_p)) + geom_histogram() + 
  facet_grid(factor~gene_set)+ ggtitle("avg logit P-values")
#It seems like for a few tests there may be some inflation- especially for the platelet function gene set
ggplot(permuted, aes(x=maxlogit_p)) + geom_histogram() + 
  facet_grid(factor~gene_set)+ ggtitle("max logit P-values")

#similarly with this one. In fact this approach seems to be the most inflated of all of them.
#How interesting. the macrogenes set