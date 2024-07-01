#Script to
unitNorm <- function(x)
{
  x/norm(x, "2")
}
unitNorms <- function(M)
{
  apply(M, 2, function(x) unitNorm(x))
}


suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("magrittr"))
suppressMessages(library("dplyr"))
suppressMessages(library(stringr))
suppressMessages(library("combinat"))
source("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/src/evaluateSimR2.R")
suppressMessages(library(NatParksPalettes))
option_list <- list(
  make_option(c("-o", "--output"), default = '', type = "character",
              help="Output dir + handle"),
  make_option(c("-s", "--sim_path"), default = '', type = "character",
              help="File containing all of the paths of simulations you wish to go into this plot."),
  make_option(c("-p", "--plot"), default = FALSE, type = "logical", action = "store_true",
              help="specify this if you want to plot."),
  make_option(c("-y", "--yaml"), default = '', type = "character",
              help="Path to the settings file."),
  make_option(c("-w", "--whitened"), default = FALSE, type = "logical",
              help="Include whitened", action = "store_true"),
  make_option(c("-n", "--scale_data"), default = FALSE, type = "logical",
              help="Scale true and predicted V/U to have unit norm.", action = "store_true")
)

#debug interactive

t = c("--output=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.3_n100gwasMF_optim_covar_high_fixed_1/factorization_results/summary",
  "--yaml=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/V1_U1_maf0.3_n100.high_covar_fixed.yml",
  "--sim_path=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.3_n100gwasMF_optim_covar_high_fixed_1/factorization_results/",  "--scale_data")
#V1_U1_maf0.01_n100.yml.gwasMF_only.yml
t = c("--output=//scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V7_U7_mafmixed_n50000.no_covar_cont_scaling//factorization_results/summary",
      "--yaml=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/V7_U7_mafmixed_n50000.no_covar_cont_scaling.yml",
      "--sim_path=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V7_U7_mafmixed_n50000.no_covar_cont_scaling//factorization_results/",  "--scale_data")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/R/sparsity_scaler.R")
#args <- parse_args(OptionParser(option_list=option_list), args = t)

args <- parse_args(OptionParser(option_list=option_list))#, args = t)

#                   args = c("--input=/Users/aomdahl/Library/CloudStorage/OneDrive-JohnsHopkins/Research_Ashton_MacBook_Pro/snp_network/scratch/cohort_overlap_exploration/babytest1"))

#Rscript src/evaluateSims.R --output simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/factorization_results/ -p -y  simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/udler_1_simple_no-corr_no-overlap.yml-s simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/factorization_results/
#Organized by directory per sim.
#corresponding script is in ldsc_all_traits/src/factMetaAssessment
#handles file, i.e. pma,pmd,flashr....
#$yaml
#[1] "./simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/udler_1_simple_no-corr_no-overlap.yml"
#args$sim.path <- "./simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/factorization_results/"
#get the names of the methods you want

yml <- read.table(args$yaml, sep = ",")
methods.run <- ((yml %>% filter(V1 == "test_methods"))$V2 %>% strsplit(.,":" ))[[1]]
true.loadings <- as.matrix(fread((yml %>% filter(V1 == "loadings"))$V2))
true.factors <-  as.matrix(fread((yml %>% filter(V1 == "factors"))$V2))
if(args$scale_data)
{
  true.loadings <- unitNorms(true.loadings)
  true.factors <- unitNorms(true.factors)
}
niter <- as.numeric((yml %>% filter(V1 == "iter"))$V2)
s <- args$sim_path
true.k <- ncol(true.factors)
#This will evaluate the simulations you desire. The assumption is that each individual simulation has its own directory, with a
#sim_loadings and sim_factors file, and then files corresponding to the predictved outcome for several different methods.
sim.performance <- NULL
f_i =1
nmethods = length(methods.run)
for(m in methods.run)
{
  r_performance <- matrix(NA, nrow = niter, ncol = 8)
  for(i in 1:niter){
    #print(i)
    #print(m)
    if(!file.exists(paste0(s, "/sim",i, ".", m, ".loadings.txt")) & i > 1)
      {
         message("WARNING: missing file: ", paste0(s, "/sim",i, ".", m, ".loadings.txt"))
         next
    }
    if(!file.exists(paste0(s, "/sim",i, ".", m, ".loadings.txt")) & i == 1)
       {
         message("WARNING: missing file: ",paste0(s, "/sim",i, ".", m, ".loadings.txt"))
         break;
    }
    #print(i)
    #print(paste0(s, "/sim",i, ".", m, ".loadings.txt"))
    pred.loadings <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".loadings.txt")))
    pred.factors <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".factors.txt")))
    if(args$scale_data)
    {
     # message("Scaling both true and loaded data for convenient comparison")
      if(!all(pred.loadings == 0))
      {
        pred.loadings <-unitNorms(pred.loadings)
        pred.factors <- unitNorms(pred.factors)
      }else
      {
        message("Unable to scale data; all elements 0'd out.")
      }

    }
    reconstruction <- fread(paste0(s, "/sim",i, ".", m, ".recon_error.txt"))

    #Get K
    ks <- ncol(pred.factors)
    sparsity.v <- matrixSparsity(pred.factors,initK = true.k, wrt.init = TRUE)
    sparsity.u <- matrixSparsity(pred.loadings,initK = true.k, wrt.init = TRUE)
    if(all(pred.loadings == 0))
    {
      ks <- 0
    }

  #maybe have the option to scale be embedded in here?
    r_performance[i,] <- c(evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.factors,unit.scale = FALSE),
                           reconstruction$Frobenius_norm[1], reconstruction$Correlation[1], ks,sparsity.v,sparsity.u, i)
    f_i = f_i+1
  }
  sim.performance <- rbind(sim.performance, data.frame("method" = m, r_performance))
}
#Now, read in whitened if they are there....
if(args$whitened)
{
  nmethods = length(methods.run) *2
  for(m in methods.run)
  {
    r_performance <- matrix(NA, nrow = niter, ncol = 5)
    for(i in 1:niter){
      if(!file.exists(paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt")) & i > 1)
      {
       message("WARNING: missing file: ", paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt"))
        next
      }
      if(!file.exists(paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt")) & i == 1)
         {
           message("WARNING: missing file: ", paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt"))
           break
      }
      pred.loadings <- fread(paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt"))
      pred.factors <- fread(paste0(s, "/sim",i, ".","whitened.", m, ".factors.txt"))
      if(args$scale_data)
      {
        message("Scaling both true and loaded data for convenient comparison")
        pred.loadings <- apply(pred.loadings, 2, function(x) x/norm(x, "2"))
        pred.factors <- apply(pred.factors, 2, function(x) x/norm(x, "2"))
      }

      reconstruction <- fread(paste0(s, "/sim",i, ".","whitened.", m, ".recon_error.txt"))
      r_performance[i,] <- c(evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.factors),
                             reconstruction$Frobenius_norm[1], reconstruction$Correlation[1], i)
      f_i = f_i+1
    }
    sim.performance <- rbind(sim.performance, data.frame("method" = paste0(m, "_whitened"), r_performance))
  }
}



colnames(sim.performance) <- c("method", 'R2_L','R2_F',"RSE", "R2_X", "K_out","sparsityV", "sparsityU", "iter")
#Now plot it, if desired
write.table(sim.performance, file = paste0(args$output, ".tabular_performance.tsv"), quote = FALSE, row.names = FALSE)
if(args$plot)
{
  #Plots of F and L
  suppressMessages(library(ggplot2))
  l.perf <- ggplot(sim.performance, aes(x = method, y = R2_L, fill = method)) + geom_boxplot() + ylab(expression(R^2)) +   ggtitle("Loadings") + theme_classic(15) +
    scale_fill_manual(values=natparks.pals("Arches", nmethods)) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank())
  f.perf <- ggplot(sim.performance, aes(x = method, y = R2_F, fill = method)) + geom_boxplot() + ylab(expression(R^2)) +
    ggtitle("Factors") + theme_classic(15) + scale_fill_manual(values=natparks.pals("Arches", nmethods))
  suppressMessages(library(cowplot))
  ggsave(plot_grid(plotlist = list(l.perf, f.perf), nrow = 2, ncol =1), filename = paste0(args$output, ".f_l_boxplot.png"))

    #plots of frobenious norm and overall R2
    frob <- ggplot(sim.performance, aes(x = method, y = RSE, fill = method)) + geom_boxplot() + ylab("Root squared error") +   ggtitle("Frobenius norm") + theme_classic(15) +
    scale_fill_manual(values=natparks.pals("Arches", (nmethods))) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank())
  r2x <- ggplot(sim.performance, aes(x = method, y = R2_X, fill = method)) + geom_boxplot() + ylab(expression(R^2)) +
    ggtitle("Overall correlation") + theme_classic(15) + scale_fill_manual(values=natparks.pals("Arches", nmethods))
  ggsave(plot_grid(plotlist = list(frob, r2x), nrow = 2, ncol =1), filename = paste0(args$output, ".recon_boxplot.png"))

  #K-plot
  ks <- ggplot(sim.performance, aes(x = method, y = K_out, fill = method)) + geom_boxplot() + ylab("Final K") +
    ggtitle("Predicted latent factors") + theme_classic(15) + scale_fill_manual(values=natparks.pals("Arches", nmethods))
  ggsave(ks, filename = paste0(args$output, ".ks_boxplot.png"))


}




