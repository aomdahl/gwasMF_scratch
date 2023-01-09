#!/usr/bin/Rscript
#renv::init("../../../custom_l1_factorization/renv_f/")
#BIC wrapper for outputs...
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("plyr")) 
suppressPackageStartupMessages(library("dplyr")) 
suppressPackageStartupMessages(library("ggplot2")) 
suppressPackageStartupMessages(library("stringr")) 
suppressPackageStartupMessages(library("penalized")) 
suppressPackageStartupMessages(library("cowplot")) 
suppressPackageStartupMessages(library("parallel")) 
suppressPackageStartupMessages(library("doParallel")) 
suppressPackageStartupMessages(library("logr")) 
suppressPackageStartupMessages(library("coop")) 
suppressPackageStartupMessages(library("data.table")) 
suppressPackageStartupMessages(library("glmnet")) 
suppressPackageStartupMessages(library("svMisc")) 
suppressPackageStartupMessages(library("nFactors")) 
#suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("optparse"))
source("/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/model_selection//bic_autofit_functions.R")
option_list <- list(
  make_option(c("--gwas_effects"), type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait"),
  make_option(c("--uncertainty"), type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait"),
  make_option(c("--lambda_gc"), type = 'character', help = "Specify the path to the genomic correction coefficients. If none provided, none used", default = ""),
  make_option(c("--trait_names"), type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the orderr in the input tables.", default = ""),
  make_option(c("--weighting_scheme"), type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE"),
  make_option(c("--init_F"), type = 'character', default = "ones_eigenvect", help = "Specify how the F matrix should be initialized. Options are [ones_eigenvect (svd(cor(z))_1), ones_plain (alls 1s), plieotropy (svd(cor(|z|))_1)]"),
  make_option(c("--init_L"), type = 'character', help = "Specify this option to start by initializing L rather than F. Options are [random, pleiotropy]. Use empty string for none (default)", default = ""),
  make_option(c("--alphas"), type = 'character', default = "", help = "Specify which alphas to do, all in quotes, separated by ',' character"),
  make_option(c("--lambdas"), type = 'character', default = "", help = "Specify which lambdas to do, all in quotes, separated by ',' character"),
  make_option(c("--scaled_sparsity"), type = "logical", action = "store_true", default = FALSE, help = "Specify this to scale the sparsity params by dataset to be between 0 and 1"),
  make_option(c("--output"), type = "character", help = "Source file location"),
  make_option(c("--calibrate_k"), type = "logical", help = "Give just an estimate of K", action = "store_true", default = FALSE),
  make_option(c("--cores"), type = "integer", help = "Number of cores", default = 1),
  make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
  make_option(c("--debug"), type = "logical", help = "if want debug run", action = "store_true", default = FALSE),
  make_option(c("--overview_plots"), type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE),
  make_option(c("-k", "--nfactors"), type = "character", help = "specify the number of factors", default = "0"),
  make_option(c("-i", "--niter"), type = "integer", help = "specify the number of iterations", default = 30),
  make_option(c("--posF"), type = "logical", default = FALSE,  help = "Specify if you want to use the smei-nonnegative setup.", action = "store_true"),
  make_option(c("--scale_n"), type = "character", default = "",  help = "Specify the path to a matrix of sample sizes if you want to scale by sqrtN as well as W"),
  make_option(c("--IRNT"), type = "logical", default = FALSE,  help = "If you want to transform the effect sizes ", action = "store_true"),
  make_option(c("--MAP_autofit"), type = "integer", default = -1,  help = "Specify if you want to the MAP autofit sparsity parameter to learn lambda, alpha"),
  make_option(c("-f", "--converged_F_change"), type="double", default=0.02,help="Change in factor matrix to call convergence"),
  make_option(c("-o", "--converged_obj_change"), type="double", default=1,help="Relative change in the objective function to call convergence"),
  make_option(c("--no_SNP_col"), type="logical", default= FALSE, action = "store_true", help="Specify this option if there is NOT first column there..."),
  make_option(c("--parameter_optimization"), type="integer", default= 1, action = "store_true", help="Specify how many iterations to run to get the cophenetic optimization, etc. "),
  make_option(c("--regression_method"), type="character", default= "penalized", help="Specify which regression method to use. Options are [penalized (Default), glmnet]"),
  make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
  make_option(c("--epsilon"), type="double", default= 1e-8, help="The convergence criteria for the L1. If exploring the space, try making this larger to speed up runtime. "),
  make_option(c("--subsample"), type="integer", default= 0, help="Specify if you wish to subsample from the full dataset, and if so by how much. For repeated runs, this will choose a different subsample each time."),
  make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud"),
  make_option(c("--auto_grid_search"), type="logical", default=FALSE, action = "store_true", help="Specify if you want to do an auatomatic grid search"),
  make_option(c("-c","--covar_matrix"), type = "character", default = "", help = "Path to a covariance matrix. Should be pre-filtered to include only significant elements"),
  make_option(c("-p","--prefix"), type = "character", default = "", help = "File prefix."),
  make_option(c("--bic_adj"), type = "integer", default = 5, help = "how many iters with BIC adjustment.")
)
args <- parse_args(OptionParser(option_list=option_list))
args <- UdlerArgs()
runFullPipeClean(args$prefix,args, gwasmfiter =args$bic_adj)
