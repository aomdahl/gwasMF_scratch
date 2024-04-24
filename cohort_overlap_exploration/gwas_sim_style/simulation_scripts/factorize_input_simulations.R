library(optparse)
#pacman::p_load(data.table, doParallel, foreach, tidyr, dplyr)
pacman::p_load(data.table, tidyr, dplyr)
option_list <- list(
  make_option(c("--input_dir"), type = 'character', help = "Specifgy the input directory with the files- assuming the newest version"),
  make_option(c("--with_covar"), type = 'logical', help = "Use flag if want to include covariance adjustment. False by default", default = FALSE, action="store_true"),
  make_option(c("-o", "--output"), type = 'character', help = "Output file path"),
  make_option(c("-s", "--start"), type = 'integer', help = "Start index", default = 1),
  make_option(c("-c", "--ncores"), type = 'numeric', help = "how many cores you wwant. if not specified, detcts (dangerous)", default=-1)
)

quiet <- function(x) { 
  sink(sink('/dev/null')) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

t <- c("--input_dir=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/large_studies_with_covar/",
       "--with_covar",
       "--output=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/large_studies_with_covar/factorization.summaries.RData")
#argv <- parse_args(OptionParser(option_list=option_list), args = t)
argv <- parse_args(OptionParser(option_list=option_list))#, args = t)
source("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/simulation_processing_tools.R")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/factorization_methods.R")

if(FALSE)
{
  #Set up parallel run
  n.cores=argv$ncores
  if(argv$ncores == -1)
  {
    n.cores <- parallel::detectCores() - 1
  }
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #how many workers are available? (optional)
  foreach::getDoParWorkers()
  
  #Stop cluster when done
  parallel::stopCluster(cl = my.cluster)
  
  
  
}
.libPaths("/home/aomdahl1/R/4.2.1")
print(.libPaths())
#assess via all the methods:
message("Starting file read in")
sim.path <- argv$input_dir
files.betas <- list.files(sim.path,pattern = "*.BETA.csv")
sim.betas <- readInGWAS(files.betas, sim.path) 
rank <- ncol(sim.betas[[1]])
message("Read in a total of ", length(sim.betas), " files containing effect size estimates.")
files.se <- list.files(sim.path,pattern = "*\\.SE.csv")
sim.se <- readInGWAS(files.se, sim.path)

#Make sure the nubmers line up
stopifnot(length(sim.se) == length(sim.betas))
stopifnot(all(sapply(files.se, 
                     function(x) gsub(x, pattern = ".SE.csv", replacement = "")) == sapply(files.betas, 
                                                                                          function(x) gsub(x, pattern = ".BETA.csv", replacement = ""))))


sim.zscores <- lapply(1:length(sim.se), function(i) as.matrix(sim.betas[[i]])/as.matrix(sim.se[[i]]))

message("Performing standard SVD on all of them:")
sim.pca.std <- pcaOnAll(sim.zscores, cols = "ALL")

message("Performing BiocSinglar PCA for model selection tools.")
sim.pca.alt <- lapply(sim.zscores, function(x) BiocSingular::runPCA(x=x, rank=rank-1, center=FALSE))
#for the eblow: this above gives the sds
#proportion var based on the PCAtools package documentation
message("Performing model selection using basic heuristics")
total.vars <- lapply(sim.zscores, function(x) sum(matrixStats::colVars(x)))
proportionvar <- lapply(1:length(sim.pca.alt), function(i) sim.pca.alt[[i]]$sdev^2/total.vars[[i]]*100) #from the PCA code
gd.point <- lapply(1:length(proportionvar), 
                   function(i) PCAtools::chooseGavishDonoho(sim.zscores[[i]], 
                                                            var.explained=proportionvar[[i]],
                                                            noise=1))
elbow.point <- lapply(proportionvar, function(x) PCAtools::findElbowPoint(x))
kaiser.point <- sapply(sim.pca.std, function(x) sum(x$d^2 > mean(x$d^2)))

if(argv$with_covar)
{
  message("Reading in covar and covar SE data")
  sim.cov <- readInGWAS(list.files(sim.path,pattern = "*\\.COV.csv"), sim.path) 
  sim.cov_se <- readInGWAS(list.files(sim.path,pattern = "*\\.COV_SE.csv"), sim.path) 
  stopifnot(length(sim.cov) == length(sim.cov_se))
  gleaner.runs.covar <- list()
}


message("Proceeding now with gleaner and flash")
flashr.runs <- list()
gleaner.runs.std <- list()

for(i in argv$start:length(sim.betas))
{
  message("On simulation number ",i )
  message(files.betas[[i]])
  suppressMessages(gleaner.runs.std[[i]] <- runSingle("GLEANER_glmnet_noCovar", as.matrix(sim.betas[[i]]), K=rank-1, se_m=as.matrix(sim.se[[i]]), covar = NULL, bic.var = "sklearn",
                                 init.mat = "V", is.sim = TRUE, 
                                 save.path = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/gleaner")) #not accounting for BIC typoe #Uhoh- we should be outputting something, otherwise the signal cutting is just too strong....
  #Works with Zou, but not with sklearn.... arg. Why why why why why why
  #Zou is a definite no-no here, I mean 8 factors is just too many.
  #Dev- zeeroes out
  suppressMessages(flashr.runs[[i]] <- flashr::flash(data=as.matrix(sim.betas[[i]])/as.matrix(sim.se[[i]]))) #uhoh!)
  
  if(argv$with_covar)
  {
    suppressMessages(gleaner.runs.covar[[i]] <- runSingle("GLEANER_glmnet", as.matrix(sim.betas[[i]]), K=rank-1, se_m=as.matrix(sim.se[[i]]), 
                                         covar = as.matrix(sim.cov[[i]]), covar_se=as.matrix(sim.cov_se[[i]]), bic.var = "sklearn",
                                       init.mat = "V", is.sim = TRUE, enforce_blocks=FALSE, shrinkWL="strimmer",
                                       save.path = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/gleaner"))#not accounting for BIC typoe #Uhoh- we should be outputting something, otherwise the signal cutting is just too strong....
    
  }
  if(i%%25 == 0)
  {
    #perform an intermediate save
    message("Intermediate update.")
    if(argv$with_covar)
    {
      save(files.betas,flashr.runs,gleaner.runs.std, gleaner.runs.covar,sim.pca.std,sim.pca.alt,gd.point,elbow.point,kaiser.point,
           file=argv$output)
    }else
    {
      save(files.betas, flashr.runs,gleaner.runs.std, sim.pca.std,sim.pca.alt,gd.point,elbow.point,kaiser.point,
           file=argv$output)
    }
  }
}

#Save the output data:
if(argv$with_covar)
{
  save(files.betas,flashr.runs,gleaner.runs.std, gleaner.runs.covar,sim.pca.std,sim.pca.alt,gd.point,elbow.point,kaiser.point,
       file=argv$output)
}else
{
  save(files.betas, flashr.runs,gleaner.runs.std, sim.pca.std,sim.pca.alt,gd.point,elbow.point,kaiser.point,
       file=argv$output)
}
#save(flashr.runs,gleaner.runs, file="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/flash_v_gleaner_N-var.RData")


