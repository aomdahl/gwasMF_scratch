########################
5/17/22
#Purpose of this script is to take an input matrix of set reference SNPS with scores of some kind and determine if they are enriched for plieotropy.
#Required input: X (original Z score data), L (factors to detect pleiotropy on), p (path to pleiotropy data), SNPs (corresponding SNP list)
########################
pacman::p_load(tidyr, dplyr, ggplot2, stringr, data.table, cowplot, magrittr, doParallel, Xmisc, logr, coop)
source("/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/helper_functions_plieo.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/pve.R")
parser <- ArgumentParser$new()
parser$add_description("Script to quickly assess the plieotropy of enriched SNPs")
parser$add_argument("--Z", type = 'character', help = "Path to file with original Z scores for testing against", default = "SEED2")
parser$add_argument("--snp_list", type = 'character', help = "List of SNPS and ", default = "SEED2")
parser$add_argument("--L_mat", type = 'character', help = "Path to L matrix. Assumes the same order as Z")
parser$add_argument("--F_mat", type = 'character', help = "Path to F matrix. Assumes the same order as Z")
parser$add_argument("--fa_reference", type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the orderr in the input tables.", 
                    default = "scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/fastAssetScores_0.1p/")
parser$add_argument("--alpha", type = 'numeric', help = "alpha threshold", 
                    default = 0.01)
parser$add_argument("--maf", type = 'character', help = "specify where the MAFs can be found. An ordered column will do.", 
                    default = "")
parser$add_argument("--output", type = 'character', help = "specify path to write ouput files", default = "")
parser$add_argument("--tab_output", type = 'logical', help = "Specify this flag to write out tabular files", action = "store_true", default = FALSE)
parser$add_argument("--debug_mode", action = "store_true", default = FALSE, help = "Specify this to run in deubg mode.", type = "logical")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
message("Please make sure the first column of input data is SNP/RSIDs.")

debug = args$debug_mode
look_path <- args$fa_reference
#look_path <- "/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/fastAssetScores_0.1p/"
if(args$Z == "SEED2" | debug)
{
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/quickLoadData.R")
  dat.gwas<- quickLoadFactorization("B_SE", "MARCC")
  trait.list <- colnames(dat.gwas$X)
  Z <- dat.gwas$X * dat.gwas$W
  snp_list <- dat.gwas$vars
  maf <- dat.gwas$MAF %>% rowMeans()
} else
{
  maf <- fread(args$maf)
}


if(debug)
{
  l.path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE.loadings.txt"
  f.path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE.factors.txt"
  args$output <- "scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE"
  
  #l.path <- "scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio_L/factorization/A2.14_L4806.343_B_SE.loadings.txt"
  #f.path <- "scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio_L/factorization/A2.14_L4806.343_B_SE.factors.txt"
  
  l.path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_gold_standard/factorization/"
  f.path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE.factors.txt"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE"
  
} else {
  l.path <- args$L_mat
  f.path <- args$F_mat
  }

message("Reading in the reference fast asset scores")
fa_outputs <- list.files(path = look_path, pattern = ".tsv")
lr <- lapply(fa_outputs, function(x) fread(paste0(look_path, x)) %>% mutate("SNP_num" = x))
fin <- do.call("rbind",lr) %>% mutate("snpid" = snp_list$ids) %>% 
  filter(cpos != 0 | cneg != 0 | pval != 0) %>% mutate("fa_fdr" = p.adjust(pval, method = "fdr"))


message("Reading in the requested data")
l.matrix <- fread(l.path) %>% filter(snps %in% unlist(snp_list$ids))
#this for plieotropy
#%>% rename("ids" = snps) %>% merge(., snp_list, by = "ids")
f.matrix <- fread(f.path)
K <- ncol(f.matrix)-1

p.values <- getPvalsRow(Z, f.matrix[,-1])
fdr.pvals <- data.frame(apply(p.values, 2, function(x) p.adjust(x, method = "fdr"))) %>% mutate("snpid" = snp_list$ids)

alpha = args$alpha
#alpha = 0.01
#FA us fast asset, nOT FACTOR ANALYSIS.
#Default setting- filter by only significant SNPs (?)
if(ncol(fdr.pvals) != (K+1))
{
    message("Error in projecting FDRs...")
}
joined.pleio <- left_join(fin, fdr.pvals, by = "snpid")  %>% 
  set_colnames(c("cpos", "cneg", "pval", "SNP_num", "snp_id","fa_fdr", paste0("FDR_", 1:K))) 
#%>% filter(fa_fdr < alpha)
joined.pleio.alpha <- joined.pleio %>% filter(fa_fdr < alpha)
message(paste0("Latent space has ", K, " factors"))
#different question- are significant L1's enriched for meta-analysis significance?

#FIRST PLOT:
#When considering SNPS significant in some trait (by Z_fa), are what is the pleiotropy score distribution of the SNPS?
factor.enrichment.joined <- NULL
t.stats <- c()
for(i in 1:K)
{
  index = 6 + i
  print(paste0("currently on factor ", colnames(joined.pleio.alpha)[index]))
  f_sig = ifelse(unlist(joined.pleio.alpha[,..index]) < alpha, "sig", "nonsig")
  df.simple <- joined.pleio.alpha %>% mutate("latent_sig" = f_sig) %>% 
    mutate("pleio_score" = cpos + cneg) %>% select(pleio_score, latent_sig, snp_id, fa_fdr, !!sym(paste0("FDR_", i)))
  
  #Test for the difference...
  r <- t.test((filter(df.simple, latent_sig == "sig"))$pleio_score, (filter(df.simple, latent_sig == "nonsig"))$pleio_score, "greater")
  t.stats <- c(t.stats, r$statistic)
  factor.enrichment.joined <- rbind(factor.enrichment.joined, df.simple %>%
                                      set_colnames(c("pleio_score", "latent_sig", "snp_id", "fa_fdr", "latent_fdr")) %>% 
                                      mutate("Factor" = paste0("F", i)))
}

factor.enrichment.joined$Factor <- factor(factor.enrichment.joined$Factor, levels = paste0("F", 1:K)[order(t.stats,decreasing = TRUE)])
ggplot(factor.enrichment.joined, aes(x = Factor, y = pleio_score, fill = latent_sig), position = position_dodge()) + labs(fill = "Latent factor\nFDR") + 
scale_fill_manual(labels = c(paste0("> ", alpha), paste0("< ", alpha)), values = c("coral", "skyblue"))+
  geom_boxplot() + theme_minimal(16) + ylab("Pleiotropy score") + ggtitle("Pleiotropy of significant SNPs by latent\nfactor (ordered by t-stat)")
ggsave(paste0(args$output, ".combined_pleiotropy.png"),width = 10)

#####
##PLOT2: Show the balance between PVE
message("Including PVE, N")
#Another option would be to order tehre by the PVE they each have. 
pma.style <- pveBySVD(as.matrix(f.matrix[,-1]),as.matrix(l.matrix[,-1]), as.matrix(Z), K = K)
flashr.style <- approxPVE_init(as.matrix(Z), as.matrix(f.matrix[,-1]),as.matrix(l.matrix[,-1]))
flashr.style2 <- approxPVE(as.matrix(f.matrix[,-1]),as.matrix(l.matrix[,-1]), as.matrix(Z), K = K)
order(pma.style, decreasing = TRUE)
order(flashr.style, decreasing = TRUE)
order(flashr.style2, decreasing = TRUE)
#perhaps a better approach- regress l directly against the ....
#testing for enrichment across all....
#How many significant SNPs per group? What is the average pleiotropy score per group?
counts.per.group <- factor.enrichment.joined %>% group_by(Factor, latent_sig) %>% summarize("n" = n(), "avg_pl" = mean(pleio_score)) %>% 
  filter(latent_sig == "sig") %>% mutate("Factor" = factor(as.character(Factor), levels =paste0("F", 1:K) )) %>% arrange(Factor)
test.df <- data.frame("pve" = pma.style, counts.per.group) %>% arrange(Factor) #The PVEs

to.plot.waste <- factor.enrichment.joined %>% left_join(.,test.df %>% select(-latent_sig), by = "Factor") %>% filter(latent_sig == "sig") %>%
  mutate("Factor" = factor(as.character(Factor), levels =paste0("F", 1:K)[order(test.df$n, decreasing = TRUE)]))
#get the order right
test.df <- test.df %>% arrange(-n)
my_xlab <- paste(test.df$Factor,"\n(N=",test.df$n,")",sep="") 
ggplot(to.plot.waste, aes(y = pleio_score, x = Factor, fill = pve)) + geom_boxplot(varwidth = TRUE) +
  scale_x_discrete(labels=my_xlab) + ylab("Pleiotropy score") + labs(fill ="PVE") + 
  theme_minimal(15) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Pleiotropy of PVE of factor SNPs (ordered by #SNPs with FDR < 0.01 per factor)")
ggsave(paste0(args$output, ".pve_sample_pleio.png"))

## Related plot: relationship between # sig SNPS and # sig traits
p.values.traits <- getPvalsCol(Z, l.matrix[,-1])
fdr.pvals.traits <- data.frame(apply(p.values.traits, 2, function(x) p.adjust(x, method = "fdr")))
sig.snps <- colSums(fdr.pvals[,-(K+1)] < alpha)
sig.traits <- colSums(fdr.pvals.traits < 0.2)

scaled.f <- apply(f.matrix[,-1], 2, function(x) x/norm(x,type = "2"))
plot.compare <- data.frame("NumSigSNPs" = sig.snps, "NumSigTraits" = sig.traits,"MinTraitFDR" = apply(fdr.pvals.traits,2,min), 
                           "Factor" = paste0("F", 1:K), "NumNonzeroTraits" = colSums(abs(scaled.f) >1e-4))
if(FALSE)
{
  ggplot(plot.compare, aes(x = NumSigSNPs, y = NumSigTraits, color = Factor)) + geom_label(aes(label = Factor)) + 
    theme_classic() + xlab(paste0("Number of significant SNPs (FDR <", alpha, ")")) + ylab("Number of signifcant traits (FDR < 0.2)")
  
  ggplot(plot.compare, aes(x = NumSigSNPs, y = MinTraitFDR, color = Factor)) + geom_label(aes(label = Factor)) + 
    theme_classic() + xlab(paste0("Number of significant SNPs (FDR <", alpha, ")")) + ylab("minimum trait FDR")

}

ggplot(plot.compare, aes(x = NumSigSNPs, y = NumNonzeroTraits, color = Factor)) + geom_label(aes(label = Factor)) + 
  theme_classic() + xlab(paste0("Number of significant SNPs (FDR <", alpha, ")")) + 
  ylab("Nonzero traits loaded in F") + theme(legend.position =  "none") 
ggsave(paste0(args$output, ".sig_snps_nonzero_traits.png"))

#########
#Significant SNP overlap between factors
overlaps <- list()
sig.snps.by.factor <- lapply(7:(K+6), function(x) joined.pleio.alpha$snp_id[fa_significant_snps[,..x] < alpha])
for(i in 1:K)
{
  overlaps[[i]] <- lapply(sig.snps.by.factor, function(x) sig.snps.by.factor[[i]] %in% x)
}
#Simple- percent sharing (does this correlate with fastAsset results?)
percent.snps.overlap <- sapply(1:K, function(i) sapply(1:K, function(j) sum(overlaps[[i]][[j]])/length(overlaps[[i]][[i]])))
#an entry in this matrix indicates what proportion of SNPS in COLUMN are also significant in ROW?
#based on my assessment here, it seems like they are generally loading on pretty different SNPs
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc//src/plot_functions.R")
plotCorrelationHeatmap(cormat =percent.snps.overlap, title = "",p.vals = NA)
colnames(percent.snps.overlap) <- paste0("F", 1:K)
rownames(percent.snps.overlap) <- paste0("F", 1:K)
order_dat <- reorder_cormat(percent.snps.overlap)
melted_cormat <- melt(order_dat$cormat)
 ggplot(data = melted_cormat, aes(x=as.factor(Var1), y=as.factor(Var2), fill=value)) + 
    geom_tile() + geom_text(aes(label=round(value * 100, digits = 1))) + theme_minimal(15) + 
   labs(fill = "% shared\nsignificant SNPs") + scale_fill_gradient(low = "white", high = "red") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Factor") + 
   ylab("Factor") + ggtitle("% overlapping SNPs between column and row")
 ggsave(paste0(args$output, ".perc_snp_overlap.png"))

#Is significance in multiple factors predictive of plieotropy?
 #as in, is a snp that is significance in multiple factors more likely to be pleiotropic?
 stopifnot(colnames(fa_significant_snps)[7] == "FDR_1")
sig.entries.per.factor <- as.data.frame(lapply(7:(K+6), function(x) fa_significant_snps[,..x] < alpha))
sig.entries.count <- rowSums(sig.entries.per.factor)
df.sig.counts <- data.frame("pleiotropy" = fa_significant_snps$cpos +  fa_significant_snps$cneg, "sig_count" = sig.entries.count)
res <- lm(pleio.score ~ sig.entries.count)
ggplot(df.sig.counts, aes(x=sig_count, y=pleiotropy)) +
  geom_jitter() + theme_classic() + 
  geom_smooth(method=lm) + xlab("Number of factors in which SNP is significant") + ylab("Pleiotropy score") + 
  annotate(geom = "text", x=max(df.sig.counts$sig_count) - 1, y = max(df.sig.counts$pleiotropy) + 1, 
           label = paste0("p = ",round(summary(res)$coef[8], digits = 36))) + 
  ggtitle("Does SNP enrichment in many factors\ncorrespond with pleiotropy?")

#######
## Which factors are predictive of plieotropy score?
message("Checking for linear enrichment now...")
cont.test <- data.frame(l.matrix) %>% rename("snpid" = snps) %>% mutate("maf" = maf) %>%
  merge(., fin, by = "snpid") %>% filter(fa_fdr < alpha) %>% mutate("pleio_score" = cpos + cneg)
joint.pleio.enrichment <- lm(cont.test$pleio_score ~ abs(as.matrix(cont.test %>% select(paste0("X", 1:K), maf))))
joint.dat <- data.frame(summary(joint.pleio.enrichment)$coefficients) %>% 
  set_colnames(c("Beta", "SE", "T", "P")) %>% mutate("Factor" = c("Intercept", paste0("F", 1:K), "MAF"))

ggplot(joint.dat[c(-1, -(K+2)),], aes(x = reorder(Factor, -T), y = T, fill = -log10(P))) + 
  geom_bar(stat = "identity") + xlab("Latent component") + ylab("T-statistic") + 
  theme_minimal(16) + ggtitle("Pleiotropy enrichment of latent factors (joint)")
ggsave(paste0(args$output, ".joint_pleiotropy.png"))

indep.tests <- lapply(1:K, function(i) summary(lm(cont.test$pleio_score ~ abs(as.matrix(cont.test %>% select(paste0("X", i), maf))))))
psandts <- data.frame(do.call("rbind", lapply(indep.tests, function(x) x$coefficients[2,c(3,4)]))) %>% 
  mutate("Factor" = paste0("F", 1:K)) %>% set_colnames(c("T", "P", "Factor"))

ggplot(psandts, aes(x=reorder(Factor, -T), y = T, fill = -log10(P))) + geom_bar(stat = "identity") + xlab("Latent component") + 
  ylab("T-statistic") + theme_minimal(16) + ggtitle("Pleiotropy enrichment of latent factors (independent)")
ggsave(paste0(args$output, ".indep_pleiotropy.png"))


###########
#Look at top loaded SNPS in l matrix- are these

#Do it by L matrix
#cumulative
scores <- NULL
t <- cont.test %>% select(pleio_score, X1) 
#first- make the column unit length
t$X1 <- t$X1/norm(t$X1, "2")
t$snp.contr.score <- (t$X1)^2 #this was a bug-- going the wrong wway.
for(a in seq(0,K, by = 0.25))
{
  f <- t %>% filter(-log10(squared.cos.score) < a) %>% select(pleio_score) %>% mutate("Score_thresh" = a)
  scores <- rbind(scores, f)
}
ggplot(scores, aes(x = Score_thresh, y = pleio_score, group = as.factor(Score_thresh))) + geom_boxplot() + theme_minimal() + xlab("SNP loading score threshold") + ylab("Pleiotropy score") + ggtitle("Pleiotropy score by SNP loading in Factor 1")
ggsave(file = paste0(args$output, ".cumulative_bins_pleio_score.png"))
#totally stratified....
for(i in 1:K)
{
  t <- cont.test %>% select(pleio_score, !!sym(paste0("X", i))) %>% set_colnames(c("pleio_score", "f"))
  t$snp.contr.score <- (t$f)^2/norm(t$f, "2")
  scores.strat <- NULL
  prev = 0
  for(a in seq(0.25,10, by = 0.25))
  {
    f <- t %>% filter(-log10(snp.contr.score) < a & -log10(snp.contr.score) > prev) %>% select(pleio_score) %>% mutate("Score_thresh" = a)
    scores.strat <- rbind(scores.strat, f)
    prev = a
  }
  ggplot(scores.strat, aes(x = Score_thresh, y = pleio_score, group = as.factor(Score_thresh))) + 
    geom_boxplot(varwidth = TRUE)+ theme_minimal() +
    geom_smooth(method = "loess", se=FALSE, color="blue", aes(group=1), linetype = "12") + xlab("Score Bin")+
    ylab("Pleiotropy score") + ggtitle(paste0("Plieotropy score across factor", i,  " score bins"))
    #ggsave(file = paste0(args$output, ".stratified_bins_factor", i, ".png")) 
}


#ACROSS FACTORS
  #TO fairly compare across factors, need to have each column scaled to unit length

scaled.l <- apply(l.matrix[,-1], 2, function(x) x/norm(x,type = "2"))

cosine.scores.l <- data.frame(t(apply(scaled.l, 1, function(x) x^2/sum(x^2)))) 
cosine.scores.l$top_factor <- apply(cosine.scores.l,1,which.max)
cosine.scores.l$snpid <- l.matrix$snps
cosine.scores.l$fdr_of_top <- sapply(1:nrow(fdr.pvals),function(i) fdr.pvals[i, cosine.scores.l$top_factor[i]])

top.pleio.snps <- fin  %>% mutate("pleio_score" = cpos + cneg) %>% select(pleio_score, snpid, fa_fdr) #%>% filter(fa_fdr < alpha)
cosine.scores.pleio.scores <- left_join(top.pleio.snps, cosine.scores.l, by = "snpid")
cosine.scores.pleio.scores$top_factor_check <- apply(cosine.scores.pleio.scores[,4:(K+3)],1, function(x) which.max(x))

stopifnot(all(cosine.scores.pleio.scores$top_factor_check == cosine.scores.pleio.scores$top_factor))

cosine.scores.pleio.scores$relative_max_contribution <- apply(cosine.scores.pleio.scores[,4:(K+3)],1, function(x) max(x))
cosine.scores.pleio.scores <- cosine.scores.pleio.scores %>% mutate("fastAsset_sig" = ifelse(fa_fdr < alpha, "sig", "nonsig"))

#Can we add information about if the SNP is significant in its own factor?

if(args$tab_output)
{
  #write the cosine scores 
  write.table(cosine.scores.pleio.scores %>% select(snpid,pleio_score,fa_fdr,fastAsset_sig,top_factor,relative_max_contribution, paste0("X", 1:K)), file = paste0(args$output, ".sig_vs_nonsig.tsv"),quote = FALSE,row.names = FALSE,sep = '\t' )
}


df.cspc <- cosine.scores.pleio.scores %>% select(pleio_score, snpid, fa_fdr, top_factor,fdr_of_top, fastAsset_sig,relative_max_contribution)
#df.cspc <- df.cspc %>% mutate("pres_score" = ifelse(fastAsset_sig == "sig", -pleio_score, pleio_score))
ggplot(df.cspc, aes(y = pleio_score, x = as.factor(top_factor), color = fastAsset_sig )) + 
  geom_jitter() + xlab("Top scoring Latent factor") + ylab("Pleiotropy score") + 
  labs(color = "fastAsset FDR < 0.01") + ggtitle("Top factor assignment of each SNP") + theme_minimal()

#It looks like the easiset thing will be just to hack a plain old miami plot...
tmp <- df.cspc %>% group_by(top_factor) %>% filter(row_number()==1) %>% 
  mutate(pleio_score = NA, fastAsset_sig = "sig", "opacity" = 0.00, "color_key" = NA) %>% filter(top_factor == 1)

#Basic miami plot
df.miami <- df.cspc %>% mutate("opacity" = 1, "color_key" = ifelse((top_factor %% 2) ==0, "skyblue", "blue")) #%>% bind_rows(., tmp)
upper <- ggplot(df.miami %>% filter(fastAsset_sig == "sig"), aes(y = pleio_score, x = as.factor(top_factor), color = color_key)) + 
  geom_jitter() + xlab("Top scoring Latent factor") + ylab(paste("FDR < ", alpha)) + 
  theme_classic(15) + scale_colour_manual(values =c("skyblue", "blue")) + 
  theme(legend.position = "none", axis.title.x = ggplot2::element_blank(), 
  plot.margin = ggplot2::margin(t = 10, b = 0, l = 10, r = 10))

lower <- ggplot(df.miami %>% filter(fastAsset_sig != "sig"), aes(y = pleio_score, x = as.factor(top_factor),color = color_key)) + 
  geom_jitter() + xlab("Top scoring Latent factor") + 
  labs(color = "fastAsset FDR < 0.01", x="") +  theme_classic(15) + scale_y_reverse() + 
  scale_colour_manual(values =c("red", "coral")) +  
  ggplot2::scale_x_discrete(position = "top") + theme(legend.position = "none", axis.title.x = ggplot2::element_blank(), 
                                                      axis.text.x=element_blank(),
                                                        plot.margin = ggplot2::margin(t = 0, b = 0, l = 10, r = 10)) + 
  ylab(paste("FDR > ", alpha))

g <- gridExtra::grid.arrange(upper, lower, nrow = 2, top = "Pleiotropy scores for top loaded factor by SNP")
ggsave(g, file = paste0(args$output, ".sig_vs_nonsig.png"))

#Plot where color is based on significance in factor
df.miami <- df.cspc %>% mutate("opacity" = ifelse(fdr_of_top < alpha, 1, 0.3), "color_key" = ifelse(fdr_of_top < alpha, "black", "grey")) #%>% 
  #bind_rows(., tmp)
upper <- ggplot(df.miami %>% filter(fastAsset_sig == "sig") %>% bind_rows(., tmp), aes(y = pleio_score, x = as.factor(top_factor), color = color_key, alpha = opacity)) + 
  geom_jitter() + xlab("Top scoring Latent factor") + ylab(paste("FDR < ", alpha)) + 
  theme_classic(15) +scale_colour_manual(values =c("black", "gray")) +  
  theme(axis.title.x = ggplot2::element_blank(), 
        plot.margin = ggplot2::margin(t = 10, b = 0, l = 10, r = 10), legend.position = "none")

lower <- ggplot(df.miami %>% filter(fastAsset_sig != "sig") %>% bind_rows(., tmp), aes(y = pleio_score, x = as.factor(top_factor),color = color_key, alpha = opacity)) + 
  geom_jitter() + xlab("Top scoring Latent factor") + 
  labs(color = "fastAsset FDR < 0.01", x="") +  theme_classic(15) + scale_y_reverse() + 
  scale_colour_manual(values =c("black", "gray"), labels = c("< 0.01", "> 0.01")) +
  ggplot2::scale_x_discrete(position = "top") + theme(legend.position = "bottom", axis.title.x = ggplot2::element_blank(), 
                                                      axis.text.x=element_blank(),
                                                      plot.margin = ggplot2::margin(t = 0, b = 0, l = 10, r = 10), ) + 
  ylab(paste("FDR > ", alpha)) +  labs(color = "Factor FDR") + guides(alpha = "none")

g <- gridExtra::grid.arrange(upper, lower, nrow = 2, top = "Pleiotropy scores for top loaded factor by SNP, \ncolored by significance within top factor")
ggsave(g, file = paste0(args$output, ".sig_vs_nonsig.with_factor_sigs.png"))

#df.cspc <- rbind(df.cspc %>% filter(fastAsset_sig == "sig"), data.frame(c(0, "na", NA, 1, "sig", NA)))

ggplot(df.cspc %>% filter(fastAsset_sig == "sig"), aes(y = pleio_score, x = as.factor(top_factor), color= relative_max_contribution)) + 
  geom_jitter() + xlab("Latent factor") + ylab("Pleiotropy score") + scale_colour_gradient2(low = "grey", high = "red") + 
  labs(color = "Relative contribution of\nfactor") + ggtitle("Relative factor weight for fastAsset FDR < 0.01 SNPs") + theme_minimal()
ggsave(paste0(args$output, ".factor_weights_FA_sig_snps.png"))

ggplot(df.cspc %>% filter(fastAsset_sig != "sig"), aes(y = pleio_score, x = as.factor(top_factor), color= relative_max_contribution)) + 
  geom_jitter() + xlab("Top scoring Latent factor") + ylab("Pleiotropy score") + scale_colour_gradient2(low = "grey", high = "red") +
  labs(color = "Relative contribution of\nfactor") + ggtitle("Relative factor weight for fastAsset FDR > 0.01 SNPs") + theme_minimal(15)
ggsave(paste0(args$output, ".factor_weight_FA_NONsig_snps.png"))

df.cspc %>% filter(fastAsset_sig != "sig") %>% group_by(top_factor) %>% summarize("n" = n())
