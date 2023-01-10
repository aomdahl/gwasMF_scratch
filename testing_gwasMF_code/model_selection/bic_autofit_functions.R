#Temp script with functions for running the BIC-based fitting approach.
#Estimator of the L1 degrees of freedom
MatrixDFV <- function(mat_in){

  sum(mat_in == 0)
}
MatrixDFU <- function(mat_in,fixed_first = FALSE)
{
  #this means we are learning U
  #Use the degrees of freedom: number of non-zero coefficients
  #sum(round(mat_in, digits = 4) != 0)
  if(fixed_first)
  {
	  sum(mat_in[,-1] != 0)
  }
  else
  {
	  sum(mat_in != 0)
  }
  #This seems like it would favor sparsity?
}


CalcMatrixBIC <- function(X,W,U,V, ev="std", weighted = FALSE, df = NULL, fixed_first = FALSE)
{
  #Calculated based on my understanding of what is given by equations 12 and 7 in Lee et al 2010
  #We assume U is known, V is our predictor here
  #11/17 UPDATE: df is calculated for U and V!
  
  n=nrow(X)
  d = ncol(X)
  if(fixed_first)
  {
	  message("Removing first factor because its fixed")
	  print(dim(V))
	  print("Above should be m x k")
	  #Regress out of X the first factor effects.
	  #Remove the first factor from the downstream steps
	  X <- X -  (U[,1] %*% t(V[,1]))
	  V <- V[,-1]
	  U <- U[,-1]
  }

  #its possible the matrices are different sizes, depending on which stage of the fitting they get passed in. Correct for this:
  #aligned <- AlignFactorMatrices(X,W,U,V); U <- aligned$U; V<- aligned$V
  k = ncol(U)
  if(!weighted)
  {
    W <- matrix(1, nrow = nrow(X), ncol = ncol(X))
  }
  if(ev =="std" )
  {
    errvar <- var(as.vector(X - U %*% t(V))) #not sure if this is right...
    message("not implemented for weighted")
  }else if(ev == "ols")
  {
    
    errvar <- AltVar(X,U)
  }else if(ev == "ols" & weighted)
  {
    errvar <- WeightedAltVar(X,W,U)
  }
  else # anuumber gets passed in.
  {
    errvar = ev
  }
  
  #This seems to be dr
  #I'm pretty sure its not.
  if(is.null(df))
  {
	  message("something went wrong,..")
    df <- MatrixDFU(V)
  }
  if(TRUE)
  {
	  ret <- norm(X*W - (U %*% t(V))*W, type = "F")^2/(n*d * errvar) + (log(n*d*k)/(n*d))*df
  }else
  {

  ret <- norm(X*W - (U %*% t(V))*W, type = "F")^2/(n*d * errvar) + (log(n*d)/(n*d))*df
  }
}

#Are you accounting for the weighting? That is important...

#Determine the variance using OLS. Basically, this is the "best" residual variance we can get.
AltVar <- function(X,U)
{
  resids <- c()
  for(col in 1:ncol(X))
  {
    fit <- lm(X[,col]~U)
    resids <- c(resids, resid(fit))
  }
  var(resids)
}


#Determine the variance using OLS. Basically, this is the "best" residual variance we can get.
WeightedAltVar <- function(X,W,U)
{
  resids <- c()
  for(col in 1:ncol(X))
  {
    w <- unlist(W[,col])
    wx <- (w*X[,col])
    wu <- w * U
    fit <- lm(wx~wu)
    resids <- c(resids, resid(fit))
  }
  
  var(resids)
}

#Helper calls
CalcVBIC <- function(X,W,U,V,fixed_first=FALSE,...)
{
  #Looking at the effect of a small change..
  #CalcMatrixBIC(X,W,U,V,df=MatrixDFV(V),...) #here, we look at sparsity, not number of non-zero coefs.
  df=MatrixDFU(V,fixed_first=fixed_first)
  CalcMatrixBIC(X,W,U,V,df=df,fixed_first = fixed_first,...)
}
CalcUBIC <- function(X,W,U,V,...)
{
  CalcMatrixBIC(t(X),t(W),V,U,df=MatrixDFU(U),...)
}


#Fit the V with the scheme:
#Iniital estimates from burn.in.sparsity and consider.parsm
FitVs <- function(X, W, initU, lambdas,option, weighted = FALSE)
{
  f.fits <- list()
  bics <- c()
  for(i in 1:length(lambdas))
  {
    l <- lambdas[[i]]
    option$lambda1 <- l
    f.fits[[i]] <- fit_F(X, W, initU, option, formerF = NULL) #%>% DropEmptyColumnsPipe(.)
  }
  #update: doing the weighted
  message("Checking the weighting")
   #For the BIC calculation, we don't  want to include dense factor 1....
  if(option$fixed_ubiq)
  {
        message("please be sure to remove F1")
  }
  if(weighted) {
    #message("Weighted way...")
    av <- WeightedAltVar(X,W,initU)
    bics <- unlist(lapply(f.fits, function(x) CalcVBIC(X,W,initU,x$V, ev=av, weighted = TRUE, fixed_first = option$fixed_ubiq)))
  }else{
    av <- AltVar(X,initU)
    bics <- unlist(lapply(f.fits, function(x) CalcVBIC(X,W,initU, x$V, ev=av, fixed_first = option$fixed_ubiq)))
  }
  
  
  return(list("fits" = f.fits, "BIC"=bics))
}

DropEmptyColumns <- function(matin)
{
  c <- colSums(matin != 0)
  if(any(c == 0))
  {
    print("dropping")
    drops <- which(c == 0)
    if(1 %in% drops){
      message("BEWARE- 1st column dropping???")
    }
    return(    matin[,-drops])
  }
  return(matin)
}

DropEmptyColumnsPipe <- function(lin)
{
  ret.lin <- lin
  ret.lin[[1]] <- DropEmptyColumns(lin[[1]])
  return(ret.lin)
  
}

#Recalculate the sparsity params for U (?)
#burn.in.sparsity <- defineSparsitySpace(X, W, option, burn.in = 5)
#consider.parsm <- SelectCoarseSparsityParams(burn.in.sparsity, 5)
#Iniital estimates from burn.in.sparsity and consider.parsm
FitUs <- function(X, W, initV, alphas,option, weighted = FALSE)
{
  l.fits <- list()
  for(i in 1:length(alphas))
  {
    a <- alphas[[i]]
    option$alpha1 <- a
    l.fits[[i]] <- fit_L(X, W, initV, option) #%>% DropEmptyColumnsPipe(.)

  }
  #TODO: recode this, so don't need the logic statement. Downstream should be able to handle it
  if(weighted) {
    av <- WeightedAltVar(t(X),t(W),initV)
    bics <- unlist(lapply(l.fits, function(x) CalcUBIC(X,W,as.matrix(x$U),initV, ev=av, weighted = TRUE)))
  }else{
    av <- AltVar(t(X),initV) #This step is quite slow.... need to speed this up somehow.
    bics <- unlist(lapply(l.fits, function(x) CalcUBIC(X,W,as.matrix(x$U),initV,ev=av)))
  }
  
  
  #return(l.fits)
  return(list("fits" = l.fits, "BIC"=bics, "resid_var" =av))
}

#From new distribution and current list, how to pick the new ones?
#Bic. list: list of BIc scores for all choices
#optimal.sparsity.param- the top parameter chosen
#new.dist: distribution of all the ne lambda parameter space
#@return a list of new sparsity points to try out.


ProposeNewSparsityParams <- function(bic.list,sparsity.params, n.points = 7, no.score = FALSE)
{
  
  if(no.score)
  {
    #then bic.list is the optimal one; generate fake scores
    fake.scores <- rep(100,length(sparsity.params))
    fake.scores[which(sparsity.params == bic.list)] <- -1
    bic.list <- fake.scores
  }
  optimal.sparsity.param <- sparsity.params[which.min(bic.list)]
  message("optimal sparsity param is ", optimal.sparsity.param)
  ordered.list <- sparsity.params[order(sparsity.params)] #order the list
  #New paradigm: always look above and below,
  #If its the smallest paramter tested
  if(min(ordered.list) == optimal.sparsity.param)
  {
    above <- ordered.list[which(ordered.list == optimal.sparsity.param) + 1]
    below <- optimal.sparsity.param - (above - optimal.sparsity.param)
    new.list <- seq(below,above,by=(above-below)/n.points)[1:(n.points-1)]
  } else if(max(ordered.list) == optimal.sparsity.param) #its the largest parameter tested
  {
    below <- ordered.list[which(ordered.list == optimal.sparsity.param) - 1]
    above <- optimal.sparsity.param + (optimal.sparsity.param-below)
    new.list <- seq(below,above,by=(above-below)/n.points)[2:(n.points)]
  }else {
    #Its bounded- our estimates should be between the one immediately above and below
    above <- ordered.list[which(ordered.list == optimal.sparsity.param) + 1]
    below <- ordered.list[which(ordered.list == optimal.sparsity.param) - 1]
    new.list <- seq(below,above,by=(above-below)/n.points)[2:(n.points-1)]
  }
  #Ensure none of them are less than 0
  if(any(new.list <= 0))
  {
    print("removing negatives...")
    rep.list <- new.list[which(new.list > 0)]
    if(any(new.list == 0))
    {
      rep.list <- c(min(rep.list)/2,rep.list)
    }
    new.list <- rep.list
  }
  c(optimal.sparsity.param, new.list)
}

quickSort <- function(tab, col = 1)
{
  tab[order(tab[,..col], decreasing = TRUE),]
}
    oneSDRule <- function(bics, params)
    {
          sd <- sd(bics)
        opt <- min(bics)
        in.range <-bics[(bics > (opt - sd)) & (bics < (opt+sd))]
        print("Optimal")
        print(opt)
        print("SD")
        print(sd)
        print(in.range)
        top.indices <- which(bics %in% in.range)
        print("all BICs")
        print(bics)
        print("tops")
        print(top.indices)
        print(params)
        optimal.l <- max(params[top.indices])
        print("Selecting params:")
        print(optimal.l)
        return(which.max(params[top.indices]))
    }

  selectOptimalInstance <- function(fit.data, bic, params, oneSD = FALSE)
  {
    if(oneSD)#just testing our the one sd rule
    {
      optimal.index <- oneSDRule(bic,params)
    } else
    {
      optimal.index <- which.min(bic)
    }
    optimal.matrix <- DropEmptyColumns(fit.data$fits[[optimal.index]][[1]])
    #SELECT NON-ZERO entries
  return(list("m" = optimal.matrix, "p" = params[[optimal.index]], "index" = optimal.index))
  }

    
    
getBICMatrices <- function(opath,option,X,W,all_ids, names, min.iter =5, max.iter = 50, burn.in.iter = 5)
{

message("Repeating number of iterations:")
#If we get columns with NA at this stage, we want to reset and drop those columns at the beginning.
  burn.in.sparsity <- DefineSparsitySpaceInit(X, W, option, burn.in = burn.in.iter) #If this finds one with NA, cut them out here, and reset K; we want to check pve here too.
  option$K <- burn.in.sparsity$new.k
  consider.params <- SelectCoarseSparsityParams(burn.in.sparsity, burn.in.iter, n.points = 15)
  
  #things to record
  rec.dat <- list("alphas"=c(), "lambdas"=c(), "bic.a" = c(), "bic.l"=c(), "obj"=c(), 
                  "v.sparsity" = c(), "u.sparsity"=c(), "iter"=c(), "sd.sparsity.u" = c(), "sd.sparsity.v" = c(),
                  "alpha.s" = c(0), "lambda.s" = c(), "Ks" = c(option$K))
  #kick things off
  lambdas <- consider.params$lambdas
  v.fits <- FitVs(X,W, burn.in.sparsity$U_burn,lambdas,option, weighted = TRUE) #also gets the bic
  optimal.iter.dat <- selectOptimalInstance(v.fits, v.fits$BIC, lambdas)
  optimal.v <- DropLowPVE(X,W,optimal.iter.dat$m) #Dropping just now...
  rec.dat$lambda.s <- c(rec.dat$lambda.s,optimal.iter.dat$p)
  #update the parameters for U based on the new V
  u.sparsity <- DefineSparsitySpace(X,W,optimal.v, "U", option)
  alphas <- SelectCoarseSparsityParamsGlobal(u.sparsity, n.points = 15) #using this because first time.
  
  #Record relevant data
  rec.dat$lambdas <- c(rec.dat$lambdas,lambdas)
  rec.dat$bic.l <- c(rec.dat$bic.l,v.fits$BIC)
  #Look at
  rec.dat$sd.sparsity.v <- c(rec.dat$sd.sparsity.v,sd(burn.in.sparsity$max_sparsity_params[[5]]$lambda))
  rec.dat$sd.sparsity.u <- c(rec.dat$sd.sparsity.u,sd(u.sparsity))
  rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
  rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v))
  #plotFactors(apply(optimal.v, 2, function(x) x / norm(x, "2")), trait_names = names, title = "best")
  NOT.CONVERGED <- TRUE; i = 1
  while(i < max.iter & NOT.CONVERGED){
    #now fit U:
	  print(i)
    u.fits <- FitUs(X, W, optimal.v, alphas,option, weighted = TRUE)
    #record relevant info- always do 
    rec.dat$alphas <- c(rec.dat$alphas, alphas);  rec.dat$bic.a <- c(rec.dat$bic.a,u.fits$BIC) 
    #Pick the best choice from here, using threshold.
    optimal.iter.dat <- selectOptimalInstance(u.fits, u.fits$BIC, alphas)
    optimal.u <- optimal.iter.dat$m
    rec.dat$alpha.s <- c(rec.dat$alpha.s,optimal.iter.dat$p)

    rec.dat$U_sparsities = c(rec.dat$U_sparsities, matrixSparsity(optimal.u, ncol(X)));

    
    #now get the new lambdas for V:
    v.sparsity <- DefineSparsitySpace(X,W,as.matrix(optimal.u),"V", option) #Is this what we want?
    rec.dat$sd.sparsity.v <- c(rec.dat$sd.sparsity.v,sd(v.sparsity))
    lambdas.new <- ProposeNewSparsityParams(v.fits$BIC, lambdas, n.points = 7)
    
    #Now fit V!
    v.fits <- FitVs(X,W, optimal.u,lambdas.new,option, weighted = TRUE)
    #Pick the best choice from here, using threshold.
    optimal.iter.dat <- selectOptimalInstance(v.fits, v.fits$BIC, lambdas.new)
    optimal.v <- DropLowPVE(X,W,optimal.iter.dat$m) #Dropping just now...
    rec.dat$lambda.s <- c(rec.dat$lambda.s,optimal.iter.dat$p)
    
    #Record new data
    rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
    rec.dat$lambdas <- c(rec.dat$lambdas, lambdas.new);  rec.dat$bic.l <- c(rec.dat$bic.l,v.fits$BIC) 

    #update the parameters for U based on the new V 
    u.sparsity <- DefineSparsitySpace(X,W,optimal.v, "U", option)
    alphas <- ProposeNewSparsityParams(u.fits$BIC, alphas,n.points = 7)
    lambdas <- lambdas.new

    rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v))
    #Check convergence
    if(i > min.iter){
      NOT.CONVERGED <- !checkConvergenceBICSearch(i, rec.dat) #returns true if convergence is reached
    }

    i = i+1
  }
  #TODO: add in drops for pve
  #This returns all the data from the last iteration
  #TODO: clean this up. This is very confusing.
  return(list("optimal.v" = optimal.v, "pve" = v.fits$fits[[optimal.iter.dat$index]]$pve,"resid.var" = u.fits$resid.var,
              "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[i], "alpha"=rec.dat$alpha.s[i], "options" = option, "K"= ncol(optimal.v)))
}
checkConvergenceBICSearch <- function(iter, record.data, conv.thresh = 1e-6)
{
  #K check:
  index = iter + 1
  queries <- c(record.data$Ks[[index]] == record.data$Ks[[index-1]],
  abs(record.data$alpha.s[[index]] - record.data$alpha.s[[index-1]]) < conv.thresh,
  abs(record.data$lambda.s[[index]] - record.data$lambda.s[[index-1]]) < conv.thresh)
  if(all(queries))
  {
    return(TRUE)
  }
  else
  {
    message("not convergned yet...")
    print(queries)
  }
  return(FALSE)
}

gwasML_ALS_Routine <- function(opath, option, X, W, optimal.v)
{
#print(W)  
reg.run <- Update_FL(X, W, option, preV = optimal.v)
  save(reg.run, file = paste0(option$out,opath, "_gwasMF_iter.Rdata" ))
print(reg.run$V)
o <- data.frame("rownames" = colnames(X), reg.run$V)
write.table(o, file =  paste0(option$out,opath, ".factors.txt"), quote= FALSE, row.names = FALSE)
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
title <- paste0("A=",reg.run$autofit_alpha[1], "Lambda=",reg.run$autofit_lambda[1]) 
plotFactors(apply(reg.run$V, 2, function(x) x/norm(x, "2")), trait_names = o$rownames, title = title)
ggsave(paste0(option$out,opath, ".factors.png"))
return(reg.run)
}

#runFullPipeClean(args$prefix,args, gwasmfiter =args$bic_adj)
runFullPipeClean <- function(opath,args, gwasmfiter =5)
{
#renv::init("../../../custom_l1_factorization/renv_f/")
  message("Current setting is to use U based sparsity each time...")
  suppressPackageStartupMessages(library("tidyr"))
  suppressPackageStartupMessages(library("plyr")) 
  suppressPackageStartupMessages(library("dplyr")) 
  suppressPackageStartupMessages(library("stringr")) 
  suppressPackageStartupMessages(library("penalized")) 
  suppressPackageStartupMessages(library("parallel")) 
  suppressPackageStartupMessages(library("doParallel")) 
  suppressPackageStartupMessages(library("logr")) 
  suppressPackageStartupMessages(library("coop")) 
  suppressPackageStartupMessages(library("data.table")) 
  suppressPackageStartupMessages(library("glmnet")) 
  suppressPackageStartupMessages(library("svMisc")) 
  suppressPackageStartupMessages(library("nFactors")) 
  suppressPackageStartupMessages(library("optparse"))
  dir ="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
  source(paste0(dir, "fit_F.R"))
  source(paste0(dir, "update_FL.R"))
  source(paste0(dir, "fit_L.R"))
  source(paste0(dir, "plot_functions.R"))
  source(paste0(dir, 'compute_obj.R'))
  source(paste0(dir, 'buildFactorMatrices.R'))
  source(paste0(dir, 'sparsity_scaler.R'))
  source(paste0(dir, 'cophenetic_calc.R'))
  source(paste0(dir, 'read_in_tools.R'))
  source(paste0(dir, 'regressionUtils.R'))
  source(paste0(dir, 'pve.R'))
  option <- readInSettings(args)
  option$swap <- FALSE
  option$alpha1 <- 1e-10
  option$lambda1 <- 1e-10
  output <- args$output
  log.path <- paste0(args$output, "gwasMF_log.", Sys.Date(), ".txt")
  lf <- log_open(log.path, show_notes = FALSE)
  options("logr.compact" = TRUE)
  options("logr.notes" = FALSE)
  
  #Read in the hyperparameters to explore
  hp <- readInParamterSpace(args)
  input.dat <- readInData(args)
  X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names
  if(option$K == 0)
  {
    message('Iniitializing X to the max')
    option$K <- ncol(X)
  }
  #Run the bic thing...
  bic.dat <- getBICMatrices(opath,option,X,W,all_ids, names)
  save(bic.dat, file = paste0(option$out,opath, "BIC_iter.Rdata" ))
  option <- bic.dat$options
  option$K <- bic.dat$K
  option$alpha1 <- bic.dat$alpha #Best from the previous. Seems to have converd..
  option$lambda1 <- bic.dat$lambda
  #option$iter <- 5
  ret <- gwasML_ALS_Routine(opath, option, X, W, bic.dat$optimal.v)
  
}


DefaultSeed2Args <- function()
{
  args <- list()
  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/B.tsv"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/SE.tsv"
  args$fixed_first <- TRUE
  args$genomic_correction <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/Lambda_gc.tsv"
  args$nfactors <- 15
  args$calibrate_k <- FALSE
  args$trait_names = ""
  args$niter <- 10
  args$alphas <- "" 
  args$lambdas <- ""
  args$autofit <- -1
  args$auto_grid_search <- TRUE
  args$cores <- 1
  args$IRNT <- FALSE
  args$weighting_scheme = "B_SE"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/model_selection/bic_autofit/"
  args$converged_obj_change <- 1
  args$scaled_sparsity <- TRUE
  args$posF <- FALSE
  args$init_F <- "ones_eigenvect"
  args$init_L <- ""
  args$epsilon <- 1e-8
  args$verbosity <- 1
  args$scale_n <- ""
  args$MAP_autofit <- -1
  args$auto_grid_search <- FALSE
  args$regression_method = "penalized"
  args$converged_obj_change <- 0.001
}

UdlerArgs <- function()
{
  args <- list()
  args$covar_matrix = ""
  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/beta_signed_matrix.tsv"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/se_matrix.tsv"
  args$fixed_first <- TRUE
  args$genomic_correction <- ""

  args$nfactors <- 57
  args$calibrate_k <- FALSE
  args$trait_names = ""
  args$niter <- 20
  args$alphas <- "" 
  args$lambdas <- ""
  args$autofit <- -1
  args$auto_grid_search <- TRUE
  args$cores <- 1
  args$IRNT <- FALSE
  args$weighting_scheme = "B_SE"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/udler_original/bic_version/"
  args$converged_obj_change <- 1
  args$scaled_sparsity <- TRUE
  args$posF <- FALSE
  args$init_F <- "ones_eigenvect"
  args$init_L <- ""
  args$epsilon <- 1e-8
  args$verbosity <- 1
  args$scale_n <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/sample_counts_matrix.tsv"
  args$MAP_autofit <- -1
  args$auto_grid_search <- FALSE
  args$regression_method = "penalized"
  args$converged_obj_change <- 0.05 #this is the percent change from one to the next.
  args$prefix <- ""
  args
}


