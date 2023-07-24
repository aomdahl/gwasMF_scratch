pacman::p_load(magrittr, tidyr, dplyr, ggplot2, data.table)
BuildYAML <- function(row)
{
  ret.list <- list()
  dir = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/"
  #factors
  ret.list[["factors"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$V, Category == "factors") %>% select(Directory, File)), collapse = "/"))
  #loadings
  ret.list[["loadings"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$U, Category == "loadings") %>% select(Directory, File)), collapse = "/"))
  #MAF
  ret.list[["maf"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$MAF, Category == "MAF") %>% select(Directory, File)), collapse = "/"))
  #K
  ret.list[["K"]] <- row$K
  #iter
  ret.list[["iter"]] <- row$NITER
  #sample size and overlap
  ret.list[["samp_overlap"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$N_o, Category == "samp_overlap") %>% select(Directory, File)), collapse = "/")) %>%
    gsub(x=., pattern = "_\\*_", replacement = paste0("_", as.character(row$N), "_"))
  
  #phenotype correlation
  ret.list[["pheno_corr"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$pheno_overlap, Category == "pheno_corr") %>% select(Directory, File)), collapse = "/")) 
  
  #test_methods
  ret.list[["test_methods"]] <- row$test_methods
  #last things (fixed)
  ret.list[["noise_scaler"]] <- 1
  ret.list[["bic_param"]] <- "sklearn"
  ret.list[["init"]] <- "V"
  if(is.na(row$HERIT))
  {
  ret.list[["herit_scaler"]] <- "continuous"

  }else
 {
	 ret.list[["herit_scaler"]] <- row$HERIT
  }


  data.frame("names" = names(ret.list), "vals" = unlist(ret.list))
}


paramfiles <- "./OneDrive - Johns Hopkins/Research_Ashton_MacBook_Pro/snp_network/simulation_to_build.csv"
conv.files <- "./OneDrive - Johns Hopkins/Research_Ashton_MacBook_Pro/snp_network/names_to_files_sims.csv"
conv.files <- "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/setting_files/path_reference.csv"
paramfiles <- "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/setting_files/v1_u1_sims.csv"
opath <- "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/25_round_sims_april23/"
param.setting <- fread(paramfiles)


param.setting <- param.setting %>% mutate("fnames" = paste0(V,"_", U, "_MAF-", MAF, "_N-",N, "_RHO-", pheno_overlap, "_No-", N_o, ".yml"))


conv <-fread(conv.files)
cats <- c("factors", "loadings", "maf", "K", "iter", "samp_overlap","pheno_corr", "test_methods","noise_scaler","bic_param", "init")

for(i in 1:nrow(param.setting))
{
	message("Currently running ", param.setting[i,]$fnames)
  tab <- BuildYAML(param.setting[i,])
  write.table(tab, file = paste0(opath, "/", param.setting[i,]$fnames), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
}


