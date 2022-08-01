## Script to simulate GWAS data from factors, loadings, and correlation structure.
#6 input files

library("optparse")
library("data.table")
option_list <- list( 
  make_option(c("-i", "--input"), default="",
              help="Path to input file handles"),
  make_option(c("-o", "--output"), default = '', type = "character", 
              help=".Output dir")
)
args <- parse_args(OptionParser(option_list=option_list))
    #args = c("--input=/Users/aomdahl/Library/CloudStorage/OneDrive-JohnsHopkins/Research_Ashton_MacBook_Pro/snp_network/scratch/cohort_overlap_exploration/babytest1"))

#read in relevant things
n_o <- as.matrix(fread(paste0(args$input, ".samp_overlap.csv")))
#n <- fread(paste0(args$input, ".samples.csv")) #on the diagonal of n_o, not needed.   
rho <- as.matrix(fread(paste0(args$input, ".pheno_corr.csv")))
f <- as.matrix(fread(paste0(args$input, ".factors.csv")))
l <- as.matrix(fread(paste0(args$input, ".loadings.csv")))
se <- as.matrix(fread(paste0(args$input, ".se.csv")))
M <- nrow(n_o)
N <- nrow(l)

#First, calculate the corelation effect matrix
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


#case of interest.... what does this look like?
#Clearly there is a strong effect we need to correct for...
#t <- svd(l %*% t(f))
#p <- svd(betas)
#image(t(t$v[,1:2]))
#image(t(p$v[,1:2]))

