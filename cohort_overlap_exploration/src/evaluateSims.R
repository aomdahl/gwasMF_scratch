#Script to 
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
              help="Path to the settings file.")
)

#debug interactive
t = c(
  "--output=/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/k2/udler1_samp-overlap_no-correlation_noise1/summary",
  "--yaml=/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/k2/udler1_samp-overlap_no-correlation_noise1/udler1_samp-overlap_no-correlation_noise1.yml",
  "--sim_path=/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/k2/udler1_samp-overlap_no-correlation_noise1/factorization_results/")
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
niter <- as.numeric((yml %>% filter(V1 == "iter"))$V2)
s <- args$sim_path

#This will evaluate the simulations you desire. The assumption is that each individual simulation has its own directory, with a 
#sim_loadings and sim_factors file, and then files corresponding to the predictved outcome for several different methods.
sim.performance <- NULL
f_i =1
for(m in methods.run)
{
  r_performance <- matrix(NA, nrow = niter, ncol = 5)
  for(i in 1:niter){
    pred.loadings <- fread(paste0(s, "/sim",i, ".", m, ".loadings.txt"))
    pred.factors <- fread(paste0(s, "/sim",i, ".", m, ".factors.txt"))
    reconstruction <- fread(paste0(s, "/sim",i, ".", m, ".recon_error.txt"))
    r_performance[i,] <- c(evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.factors),
                           reconstruction$Frobenius_norm[1], reconstruction$Correlation[1], i)
    f_i = f_i+1
  }
  sim.performance <- rbind(sim.performance, data.frame("method" = m, r_performance))
}
colnames(sim.performance) <- c("method", 'R2_L','R2_F',"RSE", "R2_X", "iter")
#Now plot it, if desired
write.table(sim.performance, file = paste0(args$output, ".tabular_performance.tsv"), quote = FALSE, row.names = FALSE)
if(args$plot)
{
  #Plots of F and L
  suppressMessages(library(ggplot2))
  l.perf <- ggplot(sim.performance, aes(x = method, y = R2_L, fill = method)) + geom_boxplot() + ylab(expression(R^2)) +   ggtitle("Loadings") + theme_classic(15) + 
    scale_fill_manual(values=natparks.pals("Arches", length(methods.run))) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank()) 
  f.perf <- ggplot(sim.performance, aes(x = method, y = R2_F, fill = method)) + geom_boxplot() + ylab(expression(R^2)) + 
    ggtitle("Factors") + theme_classic(15) + scale_fill_manual(values=natparks.pals("Arches", length(methods.run)))
  suppressMessages(library(cowplot))
  ggsave(plot_grid(plotlist = list(l.perf, f.perf), nrow = 2, ncol =1), filename = paste0(args$output, ".f_l_boxplot.png"))
  
    #plots of frobenious norm and overall R2
    frob <- ggplot(sim.performance, aes(x = method, y = RSE, fill = method)) + geom_boxplot() + ylab("Root squared error") +   ggtitle("Frobenius norm") + theme_classic(15) + 
    scale_fill_manual(values=natparks.pals("Arches", length(methods.run))) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank()) 
  r2x <- ggplot(sim.performance, aes(x = method, y = R2_X, fill = method)) + geom_boxplot() + ylab(expression(R^2)) + 
    ggtitle("Overall correlation") + theme_classic(15) + scale_fill_manual(values=natparks.pals("Arches", length(methods.run)))
  ggsave(plot_grid(plotlist = list(frob, r2x), nrow = 2, ncol =1), filename = paste0(args$output, ".recon_boxplot.png"))
  
  
  
}




