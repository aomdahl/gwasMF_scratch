########################
5/17/22
#Purpose of this script is to take an input matrix of set reference SNPS with scores of some kind and determine if they are enriched for plieotropy.
#Required input: X (original Z score data), L (factors to detect pleiotropy on), p (path to pleiotropy data), SNPs (corresponding SNP list)
########################
pacman::p_load(tidyr, dplyr, ggplot2, stringr, data.table, cowplot, magrittr, doParallel, Xmisc, logr, coop)
source("/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/helper_functions_plieo.R")
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/pve.R")
parser <- ArgumentParser$new()
parser$add_description("Script to quickly assess the plieotropy of enriched SNPs")
parser$add_argument("--Z", type = 'character', help = "Path to file with original Z scores for testing against", default = "SEED2")
parser$add_argument("--snp_list", type = 'character', help = "List of SNPS and ", default = "SEED2")
parser$add_argument("--L_mat", type = 'character', help = "Path to L matrix. Assumes the same order as Z")
parser$add_argument("--F_mat", type = 'character', help = "Path to F matrix. Assumes the same order as Z")
parser$add_argument("--fa_reference", type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the orderr in the input tables.", 
                    default = "/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/fastAssetScores_0.1p/")
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
#look_path <- "/work-zfs/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/fastAssetScores_0.1p/"
if(args$Z == "SEED2" | debug)
{
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/quickLoadData.R")
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
  l.path <- "/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE.loadings.txt"
  f.path <- "/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE.factors.txt"
  args$output <- "/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE"
  
  #l.path <- "/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio_L/factorization/A2.14_L4806.343_B_SE.loadings.txt"
  #f.path <- "/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio_L/factorization/A2.14_L4806.343_B_SE.factors.txt"
  
  l.path <- "/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_gold_standard/factorization/"
  f.path <- "/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE.factors.txt"
  args$output <- "/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/seed2_pleio/factorization/A2.14_L961.269_B_SE"
  
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
if(ncol(fdr.pvals) != K)
{
    message("Error in projecting FDRs...")
}
joined.pleio <- left_join(fin, fdr.pvals, by = "snpid")  %>% 
  set_colnames(c("cpos", "cneg", "pval", "SNP_num", "snp_id","fa_fdr", paste0("FDR_", 1:K))) %>%
  filter(fa_fdr < alpha)

message(paste0("Latent space has ", K, " factors"))
#different question- are significant L1's enriched for meta-analysis significance?

#Put all of them onto 1 plot?
factor.enrichment.joined <- NULL
t.stats <- c()
for(i in 1:K)
{
  index = 6 + i
  print(paste0("currently on factor ", colnames(joined.pleio)[index]))
  f_sig = ifelse(unlist(joined.pleio[,..index]) < alpha, "sig", "nonsig")
  df.simple <- joined.pleio %>% mutate("latent_sig" = f_sig) %>% 
    mutate("pleio_score" = cpos + cneg) %>% select(pleio_score, latent_sig, snp_id, fa_fdr, !!sym(paste0("FDR_", i)))
  
  #Test for the difference...
  r <- t.test((filter(df.simple, latent_sig == "sig"))$pleio_score, (filter(df.simple, latent_sig == "nonsig"))$pleio_score, "greater")
  t.stats <- c(t.stats, r$statistic)
  #print(ggplot(df.simple, aes(x = latent_sig, y = pleio_score, fill = latent_sig)) + geom_boxplot() +
  #  theme_minimal(15) + labs(fill = "FDR") + 
  #  scale_fill_manual(labels = c(paste0("> ", alpha), paste0("< ", alpha)), values = c("red", "blue")) + ggtitle(paste0("Pleiotropy scores in F", i, " SNPs")) + 
  #  ylab("Trait count (neg + pos)") + xlab(paste0("F", i, "SNP significance")) + 
  #  annotate(x = 2.2, y = 22, geom = "text", label = paste0("p = ",round(r$p.value, digits = 10))))
  
  factor.enrichment.joined <- rbind(factor.enrichment.joined, df.simple %>%
                                      set_colnames(c("pleio_score", "latent_sig", "snp_id", "fa_fdr", "latent_fdr")) %>% 
                                      mutate("Factor" = paste0("F", i)))
}



#HERE






#All conjoined on one plt...
factor.enrichment.joined$Factor <- factor(factor.enrichment.joined$Factor, levels = paste0("F", 1:K)[order(t.stats,decreasing = TRUE)])
ggplot(factor.enrichment.joined, aes(x = Factor, y = pleio_score, fill = latent_sig), position = position_dodge()) + labs(fill = "Latent factor\nFDR") + 
scale_fill_manual(labels = c(paste0("> ", alpha), paste0("< ", alpha)), values = c("coral", "skyblue"))+
  geom_boxplot() + theme_minimal(16) + ylab("Pleiotropy score") + ggtitle("Pleiotropy of significant SNPs by latent\nfactor (ordered by t-stat)")
ggsave(paste0(args$output, ".combined_pleiotropy.png"),width = 10)

#####
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
counts.per.group <- factor.enrichment.joined %>% group_by(Factor, latent_sig) %>% summarize("n" = n(), "avg_pl" = mean(pleio_score)) %>% 
  filter(latent_sig == "sig") %>% mutate("Factor" = factor(as.character(Factor), levels =paste0("F", 1:K) )) %>% arrange(Factor)
test.df <- data.frame("pve" = pma.style, counts.per.group) %>% arrange(Factor)
#ggplot(test.df, aes(x = n, y= pve, size = avg_pl, color = Factor)) + geom_point()
#ggplot(test.df, aes(x = Factor, y= pve, size = n, color = avg_pl)) + geom_point() + theme_minimal(16) 

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
#########

#Should I be controlling for MAF, etc. here?
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
t$squared.cos.score <- (t$X1)^2/sum((t$X1)^2)
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
  t$squared.cos.score <- (t$f)^2/sum((t$f)^2)
  scores.strat <- NULL
  prev = 0
  for(a in seq(0.25,10, by = 0.25))
  {
    f <- t %>% filter(-log10(squared.cos.score) < a & -log10(squared.cos.score) > prev) %>% select(pleio_score) %>% mutate("Score_thresh" = a)
    scores.strat <- rbind(scores.strat, f)
    prev = a
  }
  ggplot(scores.strat, aes(x = Score_thresh, y = pleio_score, group = as.factor(Score_thresh))) + 
    geom_boxplot(varwidth = TRUE)+ theme_minimal() +
    geom_smooth(method = "loess", se=FALSE, color="blue", aes(group=1), linetype = "12") + xlab("Score Bin")+
    ylab("Pleiotropy score") + ggtitle(paste0("Plieotropy score across factor", i,  " score bins"))
    ggsave(file = paste0(args$output, ".stratified_bins_factor", i, ".png")) 
}


#ACROSS FACTORS
  #TO fairly compar across factors, need to have each column scaled to unit length

scaled.l <- apply(l.matrix[,-1], 2, function(x) x/norm(x,type = "2"))
cosine.scores.l <- data.frame(t(apply(scaled.l, 1, function(x) x^2/sum(x^2)))) %>% mutate("snpid" = l.matrix$snps)


top.pleio.snps <- fin  %>% mutate("pleio_score" = cpos + cneg) %>% select(pleio_score, snpid, fa_fdr) #%>% filter(fa_fdr < alpha)
cosine.scores.pleio.scores <- left_join(top.pleio.snps, cosine.scores.l, by = "snpid")
cosine.scores.pleio.scores$top_factor <- apply(cosine.scores.pleio.scores[,3:(K+2)],1, function(x) which.max(x))
cosine.scores.pleio.scores$relative_max_contribution <- apply(cosine.scores.pleio.scores[,3:(K+2)],1, function(x) max(x))
cosine.scores.pleio.scores <- cosine.scores.pleio.scores %>% mutate("fastAsset_sig" = ifelse(fa_fdr < alpha, "sig", "nonsig"))
if(args$tab_output)
{
  #write the cosine scores 
  write.table(cosine.scores.pleio.scores %>% select(snpid,pleio_score,fa_fdr,fastAsset_sig,top_factor,relative_max_contribution, paste0("X", 1:K)), file = paste0(args$output, ".sig_vs_nonsig.tsv"),quote = FALSE,row.names = FALSE,sep = '\t' )
}


df.cspc <- cosine.scores.pleio.scores %>% select(pleio_score, snpid, fa_fdr, top_factor, fastAsset_sig,relative_max_contribution)
#df.cspc <- df.cspc %>% mutate("pres_score" = ifelse(fastAsset_sig == "sig", -pleio_score, pleio_score))
ggplot(df.cspc, aes(y = pleio_score, x = as.factor(top_factor), color = fastAsset_sig )) + 
  geom_jitter() + xlab("Top scoring Latent factor") + ylab("Pleiotropy score") + 
  labs(color = "fastAsset FDR < 0.01") + ggtitle("Top factor assignment of each SNP") + theme_minimal()

#It looks like the easiset thing will be just to hack a plain old miami plot...
tmp <- df.cspc %>% group_by(top_factor) %>% filter(row_number()==1) %>% 
  mutate(pleio_score = NA, fastAsset_sig = "sig", "opacity" = 0.00, "color_key" = NA) %>% filter(top_factor == 1)
df.miami <- df.cspc %>% mutate("opacity" = 1, "color_key" = ifelse((top_factor %% 2) ==0, "skyblue", "blue")) %>% bind_rows(., tmp)
#library(miamiplot)
#ggmiami(data = data.frame(df.miami), split_by = "fa_fdr", split_at = 0.01, 
        #p = "pleio_score", chr_colors = NULL, upper_chr_colors = c("skyblue", "blue"), lower_chr_colors =c("coral", "red"))

#no, sadly no.
#upper
upper <- ggplot(df.miami %>% filter(fastAsset_sig == "sig"), aes(y = pleio_score, x = as.factor(top_factor), alpha = opacity,color = color_key)) + 
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
