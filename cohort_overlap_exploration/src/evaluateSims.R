#Script to 
library("optparse")
library("data.table")
source("evaluateSimR2.R")
option_list <- list( 
  make_option(c("-m", "--method_handles"), default="",
              help="Path to handles associated with different methods"),
  make_option(c("-o", "--output"), default = '', type = "character", 
              help="Output dir + handle"),
  make_option(c("-s", "--sim_paths"), default = '', type = "character", 
              help="File containing all of the paths of simulations you wish to go into this plot.")
)
args <- parse_args(OptionParser(option_list=option_list))
#                   args = c("--input=/Users/aomdahl/Library/CloudStorage/OneDrive-JohnsHopkins/Research_Ashton_MacBook_Pro/snp_network/scratch/cohort_overlap_exploration/babytest1"))


#Organized by directory per sim.
#corresponding script is in ldsc_all_traits/src/factMetaAssessment
#handles file, i.e. pma,pmd,flashr....
#args$method_handles <- "/Users/aomdahl/Library/CloudStorage/OneDrive-JohnsHopkins/Research_Ashton_MacBook_Pro/snp_network/scratch/cohort_overlap_exploration/method.handles.tmp"
#args$true_factors <- "/Users/aomdahl/Library/CloudStorage/OneDrive-JohnsHopkins/Research_Ashton_MacBook_Pro/snp_network/scratch/cohort_overlap_exploration/babytest1/babytest1"
#args$simdir <- "/Users/aomdahl/Library/CloudStorage/OneDrive-JohnsHopkins/Research_Ashton_MacBook_Pro/snp_network/scratch/cohort_overlap_exploration/babytest1"

file.handles <- scan(args$method_handles, what = character())
sim.paths <- scan(args$sim_paths, what = character())
#This will evaluate the simulations you desire. The assumption is that each individual simulation has its own directory, with a 
#sim_loadings and sim_factors file, and then files corresponding to the predictved outcome for several different methods.
sim.performance <- NULL
for(s in sim.paths)
{
  true.loadings <- as.matrix(fread(paste0(s, "sim_loadings.csv")))
  true.factors <-  as.matrix(fread(paste0(s, "sim_factors.csv")))
  r_performance <- matrix(NA, nrow = length(file.handles), ncol = 2)
  f_i =1
  for(f in file.handles)
  {
    pred.loadings <- fread(paste0(s, "/", f, ".loadings.txt"))
    pred.factors <- fread(paste0(s, "/", f, ".loadings.txt"))
    r_performance[i,] <- evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.factors)
    f_i = f_i+1
  }
  sim.performance <- rbind(sim.performance, data.frame("method" = f, r_performance, "sim" = s))
}

#Now plot it, if desired
write.table(sim.performance, file = paste0(args$output, ".tabular_performance.tsv"), quote = FALSE, row.names = FALSE)
