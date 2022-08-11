## Script to simulate GWAS data from factors, loadings, and correlation structure.
#6 input files
varReport <- function(t, m, n)
{
  return(c(var(t), var(m)/var(t), var(n)/var(t)))
}

library("optparse")
library("data.table")
library("ggplot2")
library("magrittr")
option_list <- list( 
  make_option(c("-i", "--input"), default="",
              help="Path to input file handles"),
  make_option(c("-s", "--seed"), default=22,
              help="Set random seed, default is 22"),
  make_option(c("-o", "--output"), default = '', type = "character", 
              help="Output dir")
)
#args <- parse_args(OptionParser(option_list=option_list), args = t)
args <- parse_args(OptionParser(option_list=option_list))

t = c("--input=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/k2/udler1_realistic-high_all-correlated-weak_noise1/udler1_realistic-high_all-correlated-weak_noise1.yml")

#read in relevant things
#Maybe a better way would be to specify a parameter file, that has all of this information in it.
#This way requires duplicates of each file in every simulation directory, which isn't what we want, is it.
#Specify a setup file.
yml <- read.table(args$input,header = FALSE,sep = ",") %>% set_colnames(c("n", "p"))
n_o <- as.matrix(fread(unlist(yml[which(yml$n == "samp_overlap"),2])))
rho <- as.matrix(fread(unlist(yml[which(yml$n == "pheno_corr"),2])))
f <- as.matrix(fread(unlist(yml[which(yml$n == "factors"),2])))
l <-  as.matrix(fread(unlist(yml[which(yml$n == "loadings"),2])))
noise.scaler = 1
if("noise_scaler" %in% yml$n)
{
	noise.scaler = as.numeric(yml[which(yml$n == "noise_scaler"),2])
}
message(paste0("Noise scaler is ", noise.scaler))
maf.mat <- as.matrix(fread(unlist(yml[which(yml$n == "maf"),2])))
M <- nrow(n_o)
N <- nrow(l)

#Set up the necessary data for SE estimation
n.mat <- matrix(do.call("rbind", lapply(1:N, function(x) diag(n_o))), nrow = N, ncol = M)
var.maf <- 2*(maf.mat)*(1-maf.mat)
S <- 1/sqrt(var.maf *n.mat)
#maf.mat <- matrix(rep(maf$maf[1:500],10), nrow= 500, ncol = 10)

#Quick sanity check: do our estimated S correspond with the actual S?
if("se" %in% yml$n)
{
  se <- as.matrix(fread(unlist(yml[which(yml$n == "se"),2])))
  s.p <- data.frame("est" = S[,1], "t"=se[,1])
  ggplot(s.p, aes(x = est, y = t)) + geom_point() + theme_classic(15) + ylab("Actual SE")+ xlab("Sim")
  suppressMessages(ggsave(filename = paste0(args$output, ".se_sim_check.png")))
}


#First, calculate the correlation effect matrix
#This should be symmetric...
#if its not, that's an error on the part of
#C is m x m
C <- matrix(NA, nrow = M, ncol = M)
for(i in 1:nrow(n_o))
{
  for(j in 1:ncol(n_o))
  {
    C[i,j] <- (n_o[i,j]/(sqrt(n_o[i,i]) * sqrt(n_o[j,j]))) * rho[i,j]
  }
}
stopifnot(isSymmetric(C))
#NOTE: could also implement in terms of the matrix normal, a single line. That would work too, but I 
#think (?) would be the same

#Create an image with the correlation structure:
source("/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
o <- plotCorrelationHeatmap(C, typin = "None", title = "Generated matrix of spurious covariance structure.")
suppressMessages(ggsave(plot = o, filename = paste0(args$output, ".spurious_covariance.png")))
write.table(x = C, file = paste0(args$output, ".c_matrix.txt"), quote = FALSE, row.names = FALSE)

#Next, generate all the effect sizes and "true" sigmas.
suppressMessages(library(MASS))
betas <- matrix(NA, nrow = N, ncol = M)
out.ses <- matrix(NA, nrow = N, ncol = M)
var.dat <- matrix(NA, nrow = N, ncol = 3)

if(FALSE) #calibrating the noise scaler for this run..
{
  for(i in 1:100)
  {
    s = 10
    mu <- f %*% l[s,]
    #S <- diag(se[s,]) %*% solve(C)
    #sigma <- sqrt(diag(se[s,])) %*% C %*% sqrt(diag(se[s,]))
    #as of 8/9
    s.sqrt <- (diag(S[s,])) #dropped the square root, the scaling was off and we want variance, not standard deviation
    sigma <- s.sqrt %*% C %*% s.sqrt
    #2 ways this could be done.....
    #betas[s,] <- mvrnorm(n = 1, mu, 5 * sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    #Or
    noise <- i*mvrnorm(n = 1, rep(0, length(mu)), 5 * sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    betas[s,] <- mu + noise
    var.dat[i,] <- varReport(betas[s,], mu, noise)
  }
  plot(1:100, var.dat[1:100,3], main = "Prop. variance in noise")
  plot(1:100, var.dat[1:100,2], main = "Prop. Variance in data")
}
for(s in 1:N)
{
  mu <- f %*% l[s,]
  #S <- diag(se[s,]) %*% solve(C)
  #sigma <- sqrt(diag(se[s,])) %*% C %*% sqrt(diag(se[s,]))
  #as of 8/9
  s.sqrt <- (diag(S[s,])) #dropped the square root, the scaling was off and we want variance, not standard deviation
  sigma <- s.sqrt %*% C %*% s.sqrt
  #2 ways this could be done.....
  #betas[s,] <- mvrnorm(n = 1, mu, 5 * sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  #Or
  noise <- noise.scaler*mvrnorm(n = 1, rep(0, length(mu)), 5 * sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  betas[s,] <- mu + noise
  #varReport(betas[s,], mu, noise)
  
  var.dat[s,] <- varReport(betas[s,], mu, noise)
  
  out.ses[s, ] <- diag(sigma) #variance of individual betas; context is now missing. Seems a bit.. odd.
  #With this mind, we could just be using the sigma hats directly, don't need to do this little dance.
  #waste of fetching time, that. Not really though, I learned some new things
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
out.beta <- data.frame("SNP" = paste0("rs", 1:N), round(betas, digits = 7)) %>% set_colnames(c("SNP", paste0("T", 1:M)))
#TODOIPDATE
out.se <- data.frame("SNP" = paste0("rs", 1:N), round(out.ses, digits = 7)) %>% set_colnames(c("SNP", paste0("T", 1:M)))


#Write out output:
write.table(x = out.beta, file = paste0(args$output, ".effect_sizes.txt"), quote = FALSE, row.names = FALSE)
write.table(x = out.se, file = paste0(args$output, ".std_error.txt"), quote = FALSE, row.names = FALSE)

full_cover <- S %*% C %*% t(S)
write.table(x = C, file = paste0(args$output, ".c_matrix.txt"), quote = FALSE, row.names = FALSE)
sink(paste0(args$output, ".variance_notes.txt"))
cat(paste0("On average, ", round(mean(var.dat[,2])*100, digits =2), "% of sample variation from true value.\n"))
cat(paste0("On average, ", round(mean(var.dat[,3])*100,digits = 2), "% of sample variation from noise."))
sink()

message(paste0("On average, ", round(mean(var.dat[,2])*100, digits =2), "% of sample variation from true value."))
message(paste0("On average, ", round(mean(var.dat[,3])*100,digits = 2), "% of sample variation from noise."))

#need a good way to look at SEs....

#The 2nd PC correlates strongly with sample size
#the 1st PC correlates strongly with MAF?


