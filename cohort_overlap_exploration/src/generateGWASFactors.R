## Script to simulate GWAS data from factors, loadings, and correlation structure.
#6 input files

library("optparse")
library("data.table")
library("ggplot2")
library("magrittr")
option_list <- list( 
  make_option(c("-i", "--inputs"), default="",
              help="Path to input file handles"),
  make_option(c("-s", "--seed"), default=22,
              help="Set random seed, default is 22"),
  make_option(c("-o", "--output"), default = '', type = "character", 
              help=".Output dir")
)
args <- parse_args(OptionParser(option_list=option_list))
set.seed(args$seed)
#read in relevant things
#Maybe a better way would be to specify a parameter file, that has all of this information in it.
#This way requires duplicates of each file in every simulation directory, which isn't what we want, is it.
#Specify a setup file.
yml <- fread(args$inputs,header = FALSE) %>% set_colnames(c("n", "p"))
print(yml[which(yml$n == "samp_overlap"),2])
n_o <- as.matrix(fread(unlist(yml[which(yml$n == "samp_overlap"),2])))

rho <- as.matrix(fread(unlist(yml[which(yml$n == "pheno_corr"),2])))
f <- as.matrix(fread(unlist(yml[which(yml$n == "factors"),2])))
l <-  as.matrix(fread(unlist(yml[which(yml$n == "loadings"),2])))
se <- as.matrix(fread(unlist(yml[which(yml$n == "se"),2])))
M <- nrow(n_o)
N <- nrow(l)

#First, calculate the correlation effect matrix
#C is m x m
C <- matrix(NA, nrow = M, ncol = M)
for(i in 1:nrow(n_o))
{
  for(j in 1:ncol(n_o))
  {
    C[i,j] <- (n_o[i,j]/(sqrt(n_o[i,i]) * sqrt(n_o[j,j]))) * rho[i,j]
  }
}
#NOTE: could also implement in terms of the matrix normal, a single line. That would work too, but I 
#think (?) would be the same

#Next, generate all the effect sizes
library(MASS)
betas <- matrix(NA, nrow = N, ncol = M)
for(s in 1:N)
{
  mu <- f %*% l[s,]
  
  sigma <- sqrt(diag(se[s,])) %*% C %*% sqrt(diag(se[s,]))
  betas[s,] <- mvrnorm(n = 1, mu, sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
}


#Scale so z-scores are z-scores? Not sure if this is necessary, but might be worth including.
#If my initial estimates of u, v are centered at 0, this shouldn't be necessary
if(FALSE)
{
  hist(out.beta[,3]/out.se[,3])
  z_tilde <- betas/se
  sigma_z <- sqrt(do.call("rbind", lapply(1:N, function(i) apply(z_tilde, 2, var))))
  
  meanz <- do.call("rbind", lapply(1:N, function(i) apply(z_tilde, 2,mean)))
  
  #beta_tilde <- (betas / sigma_z) - (meanz * se/sigma_z)
  beta_tilde <- betas / sqrt(do.call("rbind", lapply(1:N, function(i) apply(betas, 2, var))))
  #I need oto think about this a little more., Not sure how well things will work if z scores aren't N(0,1).
  #Maybe a non issue, idk.
  
}

library(magrittr)
out.beta <- data.frame("SNP" = paste0("rs", 1:N), round(betas, digits = 5)) %>% set_colnames(c("SNP", paste0("T", 1:M)))
out.se <- data.frame("SNP" = paste0("rs", 1:N), round(se, digits = 5)) %>% set_colnames(c("SNP", paste0("T", 1:M)))


#Write out output:
write.table(x = out.beta, file = paste0(args$output, ".effect_sizes.txt"), quote = FALSE, row.names = FALSE)
write.table(x = out.se, file = paste0(args$output, ".std_error.txt"), quote = FALSE, row.names = FALSE)

#Create an image with the correlation structure:
source("/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
o <- plotCorrelationHeatmap(C, typin = "None", title = "Generated matrix of spurious covariance structure.")
ggsave(plot = o, filename = paste0(args$output, ".spurious_covariance.png"))


