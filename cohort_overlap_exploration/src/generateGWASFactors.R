## Script to simulate GWAS data from factors, loadings, and correlation structure.
#6 input files

#Returns the variance of the betas (overal variance)
#variance of true data of the total
#variance of the noise given the total.
#t is betas
#m is mean
#n is noise
varReport <- function(effects, true.signal, noise)
{
  t <- unlist(as.vector(effects))
  m <- unlist(as.vector(true.signal))
  n <- unlist(as.vector(noise))
  stopifnot(length(t) == length(n))
  #total variance, fraction from true signal, fraction from noise, var true /var noise, #var mu over sum var #var noise over sum var
  return(c(var(t), var(m)/var(t), var(n)/var(t), var(m)/var(n), var(m)/(var(n) +  var(m)),var(n)/(var(n) +  var(m)) ))
}

library("optparse")
library("data.table")
library("ggplot2")
library("magrittr")
option_list <- list(
  make_option(c("-i", "--input"), default="",
              help="Path to input yaml file"),
  make_option(c("-s", "--seed"), default=22,
              help="Set random seed, default is 22"),
  make_option(c("-o", "--output"), default = '', type = "character",
              help="Output dir")
)
#t = c("--input=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/k4/noise2/udler2_realistic-high-1_r-75_noise2/udler2_realistic-high-1_r-75_noise2.yml")
t = c("--input=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/updated_variance_scaling/V6_U6_mafmixed_n100000.high_covar_1block_cont_scaling.yml",
      "--output=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs//V6_U6_mafmixed_n100000.high_covar_1block_cont_scaling/fake")

#args <- parse_args(OptionParser(option_list=option_list), args = t)
args <- parse_args(OptionParser(option_list=option_list))
#set.seed(args$seed)

#read in relevant things
#Maybe a better way would be to specify a parameter file, that has all of this information in it.
#This way requires duplicates of each file in every simulation directory, which isn't what we want, is it.
#Specify a setup file.
yml <- read.table(args$input,header = FALSE,sep = ",") %>% set_colnames(c("n", "p"))
if("sim_ref" %in% yml$n)
{
  message("Simulation path data already provided; Will skip simulation generation")
  quit()
}


n_o <- as.matrix(fread(unlist(yml[which(yml$n == "samp_overlap"),2])))

new.seed=args$seed + as.numeric(n_o[1,1])
message("New sim version: seed set by setting + n: ", new.seed)
set.seed(new.seed)
rho <- as.matrix(fread(unlist(yml[which(yml$n == "pheno_corr"),2]))) #Correlation between phenotypes
f <- as.matrix(fread(unlist(yml[which(yml$n == "factors"),2])))
l <-  as.matrix(fread(unlist(yml[which(yml$n == "loadings"),2]))) #gut says this should be more dense, huh.
noise.scaler = 1
if("noise_scaler" %in% yml$n)
{
	noise.scaler = as.numeric(yml[which(yml$n == "noise_scaler"),2])
}
message(paste0("Noise scaler is ", noise.scaler))
maf.mat <- as.matrix(fread(unlist(yml[which(yml$n == "maf"),2])))
M <- nrow(n_o)
N <- nrow(l)

if("herit_scaler" %in% yml$n)
{
  #Has to correspond to the number of features
  ntraits = nrow(f)
  message("Scaling X according to empirical trait heritabilities")
  herit.settings <- unlist(yml[which(yml$n == "herit_scaler"),2])
  
  #Effect size distribution heritabilities come from 2018 text:
  #"Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits"
  #Supplementary Table 4
  dis.traits <- c("Alzheimers","Asthma","Coronary artery disease", "Type 2 diabetes","Crohns disease", "Inflammatory bowel disease",
                      "Ulcerative colitis","Rheumatoid arthritis")
  dis.herit <- c(173.7,400.6,134.6,198.4,380.5,214.3,301.3,174.6)*10^-5
  
  continuous_traits <- c("Age at menarche","BMI","Height","Hip circumference","Waist circumference","Waist-to-hip ratio","HDL cholesterol",
                         "LDL cholesterol","Total cholesterol","Triglycerides","Child birth weight","Childhood obesity","Neuroticism")
  
  cont.herit <- c(12.4,18.0,14.6,8.0,8.1,6.5,15.2,21.8,22.7,23.1,9.1,95.6,0.9)*10^-5
  
  #quick.test

  #choose the one based on the herit.settings
  herit.factors <- switch(herit.settings,
         "disease" = rep(dis.herit,5)[1:ntraits],
         "mixed" = c(rep(dis.herit,5)[1:floor(ntraits/2)],rep(cont.herit,5)[1:floor(ntraits/2)]),
         "continuous" = rep(cont.herit,5)[1:ntraits],
         "ukbb_real" = rep(variance.scale,5)[1:ntraits])
  scaling.d <- diag(herit.factors)
  #need to multiply each column of mu by  sqrt(scaling.d[i]/var(i))
  
  #TODO- implement how this D gets applied!
  
  
}


#Set up the necessary data for SE estimation
n.mat <- matrix(do.call("rbind", lapply(1:N, function(x) diag(n_o))), nrow = N, ncol = M)
var.maf <- 2*(maf.mat)*(1-maf.mat)
S <- 1/sqrt(var.maf *n.mat) #This is the standard error
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
library(matrixcalc)
if(!is.positive.definite(C))
{
  message("C is not PD, adjusting now...")
  library(corpcor)
  C.updated <- make.positive.definite(C, tol = 1e-3)
  #percent change
  prop.change <- norm(C- C.updated, type = "F")/norm(C, type = "F")
  C <- C.updated
}
#NOTE: could also implement in terms of the matrix normal, a single line. That would work too, but I
#think (?) would be the same

#Create an image with the correlation structure:
source("/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
rownames(C) <- paste0("T", 1:nrow(C)); colnames(C) <- paste0("T", 1:ncol(C))

o <- plotCorrelationHeatmap(C, typin = "None", title = "Generated matrix of spurious covariance structure.", show.nums = TRUE)
suppressMessages(ggsave(plot = o, filename = paste0(args$output, ".spurious_covariance.png")))
write.table(x = C, file = paste0(args$output, ".c_matrix.txt"), quote = FALSE, row.names = FALSE)

#Next, generate all the effect sizes and "true" sigmas.
suppressMessages(library(MASS))
betas <- matrix(NA, nrow = N, ncol = M)
out.ses <- matrix(NA, nrow = N, ncol = M)
var.dat <- matrix(NA, nrow = N, ncol = 6)

#reset the seed- not sure if the above matrix functions require random seed, but just in case:
global.noise <- c()
global.signal <- c()
set.seed(args$seed)
all.noise <- NULL
mu.tot <- l %*% t(f)
#Need to reorganize how i do the sims. turd.
for(s in 1:N) #N is number of SNPs
{
  mu <- f %*% l[s,]
  #Variance on beta is standard error squared.
  se.i <- (diag(S[s,]))
  sigma <- se.i %*% C %*% se.i
  #2 ways this could be done.....
  #betas[s,] <- mvrnorm(n = 1, mu, 5 * sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  #Or
  #What is the tolerance argument?

  #noise <- noise.scaler*mvrnorm(n = 100, as.vector(rep(0, length(mu))), sigma, empirical = TRUE) #tol = 1e-3,
  noise <- (noise.scaler * mvtnorm::rmvnorm(1, sigma = sigma))
  all.noise <- rbind(all.noise, noise)
  global.noise <- c(global.noise, unlist(as.vector(noise)))
  global.signal <- c(global.signal, unlist(as.vector(mu)))
  #When compared to a rnorm version when indpendent samples, this is the same. so that's fine.
  betas[s,] <- mu + t(noise)

  #Write this out somehwere useful.
  var.dat[s,] <- varReport(betas[s,], mu, noise)

  out.ses[s, ] <- diag(sigma) #variance of individual betas; context is now missing. Seems a bit.. odd.
  #With this mind, we could just be using the sigma hats directly, don't need to do this little dance.
  #waste of fetching time, that. Not really though, I learned some new things
}

#Alternative version, where we add all the noise at the end
betas.alt <- mu.tot + all.noise
#verify
#stopifnot(betas == betas.alt)

#Now scale mu to match our distributional assumptions
if("herit_scaler" %in% yml$n)
{
  #May updates- only scaling the active SNPs to match the desired distribution:
  #get the scaling factors
  mu.tot.var <- apply(mu.tot, 2, function(x) var(x[x!=0]))
  #mu.tot.var <- apply(mu.tot, 2, var)
  #scaling.mat <- diag(sqrt(herit.factors/mu.tot.var))
  #mu.scaled <- (mu.tot %*% scaling.mat) 

  mu.scaled <- mu.tot
  
  for(i in 1:length(mu.tot.var)) {mu.scaled[(mu.scaled[,i] !=0),i] <- mu.scaled[(mu.scaled[,i] !=0),i] * sqrt(herit.factors[i]/mu.tot.var[i])} #scaling just the variants that are active
  #sanity check
  new.vars <- apply(mu.scaled, 2, function(x) var(x[x != 0]))
  stopifnot(max(new.vars - herit.factors) < 1e-10)
  betas = mu.scaled + all.noise
  #Some plots if doing this manually
  
  pb <- cor(betas); rownames(pb) = paste0("T", 1:10); colnames(pb) = paste0("T", 1:10)
  bcor <- plotCorrelationHeatmap(pb,typin = "None",title = "Correlation structure of scaled beta hats")
  suppressMessages(ggsave(plot = bcor, filename = paste0(args$output, ".effect_size_estimate_correlation.png")))
  if(FALSE){
    plot(herit.factors, new.vars,pch=19, xlab = "True var (GWAS)", ylab = "Scaled var"); abline(a=0, b= 1, col = "blue")
    pba <- cor(betas.alt); rownames(pba) = paste0("T", 1:10); colnames(pba) = paste0("T", 1:10)
    plotCorrelationHeatmap(pba,typin = "None",title = "Correlation structure of original beta hats")
  }
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
out.var <- data.frame("SNP" = paste0("rs", 1:N), round(var.dat, digits = 7)) %>%
  set_colnames(c("SNP", "var_beta", "var_mu:var_beta","var_noise:var_beta", "var_mu:var_noise", "var_m:var_m+var_n", "var_n:var_m+var_n"))

##Added in to run factorGo- write out Z and N:
out.z <- data.frame("SNP" = paste0("rs", 1:N), round(betas/out.ses, digits = 7)) %>% set_colnames(c("SNP", paste0("T", 1:M)))
out.n <- data.frame("N"=colMeans(n.mat))
#total variance, fraction from true signal, fraction from noise, var true /var noise, #var mu over sum var #var noise over sum vara
#Write out output:
write.table(x = out.beta, file = paste0(args$output, ".effect_sizes.txt"), quote = FALSE, row.names = FALSE)
write.table(x = out.se, file = paste0(args$output, ".std_error.txt"), quote = FALSE, row.names = FALSE)
write.table(x=var.dat, file = paste0(args$output, ".variance_data.txt"), quote = FALSE, row.names = FALSE)
full_cover <- S %*% C %*% t(S)
write.table(x = C, file = paste0(args$output, ".c_matrix.txt"), quote = FALSE, row.names = FALSE)
#write out empirical covar
write.table(x = cor(all.noise), file = paste0(args$output, ".empirical_covar.txt"), quote = FALSE, row.names = FALSE)
global.var.report <- varReport(betas, global.signal, global.noise)

#Write out for factorGo
write.table(x=out.z, file = paste0(args$output, ".z.txt"), quote = FALSE, row.names = FALSE,sep = "\t")
write.table(x=out.n, file = paste0(args$output, ".N.txt"), quote = FALSE, row.names = FALSE)


sink(paste0(args$output, ".variance_notes.txt"))
cat(paste0("Globally, ", round(mean(global.var.report[5])*100, digits =2), "% of sample variation from true value.\n"))
cat(paste0("Globally, ", round(mean(global.var.report[6])*100,digits = 2), "% of sample variation from noise."))
cat(paste0("Globally, signal to noise ratio is:", round(mean(global.var.report[4]),digits = 4)))
cat(paste0("On average, ", round(mean(var.dat[,5])*100, digits =2), "% of sample variation from true value.\n"))
cat(paste0("On average, ", round(mean(var.dat[,6])*100,digits = 2), "% of sample variation from noise."))
cat(paste0("On average, signal to noise ratio is:", round(mean(var.dat[,4]),digits = 4)))
sink()

message(paste0("Globally, ", round(mean(global.var.report[5])*100, digits =2), "% of sample variation from true value.\n"))
message(paste0("Globally, ", round(mean(global.var.report[6])*100,digits = 2), "% of sample variation from noise."))
message(paste0("On average, ", round(mean(var.dat[,5])*100, digits =2), "% of sample variation from true value."))
message(paste0("On average, ", round(mean(var.dat[,6])*100,digits = 2), "% of sample variation from noise."))
message(paste0("On average, signal to noise ratio is:", round(mean(var.dat[,4]),digits = 4)))
#need a good way to look at SEs....

#The 2nd PC correlates strongly with sample size
#the 1st PC correlates strongly with MAF?


