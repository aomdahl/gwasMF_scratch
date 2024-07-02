#Condensed code extending functionality explored in `comparison_methods.Rmd`
#The purpose of this script is to compare 2 GWAS factorizations for similarity, using standard 
#ML-style accuracy scores.


#setup
source("~/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc//src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/src/evaluateSimR2.R")

pacman::p_load(ggplot2, magrittr, tidyr)

########Helper functions


#'Get the correlation via procustes, standard CARET, PR metrics of 2 matrices
#'
#' @param A the "reference matrix"
#' @param B the 
#'
#' @return list object containing matrices with their columns ordered, mapped via procstes analysis, CARET package reports, PR curve metrics 
#' @export

#For testing: 
#A <- matrix(rnorm(1000), nrow = 100)
#B <- matrix(rnorm(1000), nrow = 100)
compareModelMatricesComprehensive <- function(A,B, corr.type = "pearson", full.procrust=TRUE)
{
  library(PRROC)
  SWAP = FALSE #refers to sizing order to get to the same size.
  lead <- as.matrix(A); scnd <- as.matrix(B)
  if(ncol(B) > ncol(A))
  {
    SWAP = TRUE
    lead <- as.matrix(B)
    scnd <- as.matrix(A)
  }
  #Get the match and the greedy correlation corresponding to it
  aligned.corr <- pseudoProcrustes(lead, scnd,corr.type = corr.type) #TODO:check this is reasonably optimal.
  lead <- aligned.corr$A
  scnd <- aligned.corr$B
  #alternative: kendall:
  aligned.corr.k <- pseudoProcrustes(lead, scnd,corr.type = "kendall") #TODO:check this is reasonably optimal. #this is much harder to calculate for U, weird still slow
  #Do a classical procrutes analysis and get the correlation
  cor.by.procrust <- NA; rss.by.procrust <- NA
  if(full.procrust)
  {
    procestes.corr <- fullProcrustes(lead, scnd) #This is too slow.
    cor.by.procrust <- stackAndAssessCorr(procestes.corr$A.normalized, procestes.corr$B.transformed)
    rss.by.procrust <- procestes.corr$RSS
  }

  #message("### Using pseudoProcrustes to calculate correlation for speed, change this in final run.")
  pseudo.procestes.corr <- aligned.corr$greedy.match$corr
  
  #do a simple BINARY precision variance assessment
  prd <-getBinaryPRCurve(lead, scnd)

  #Do a multi-class assessment, looking at +/0/-
  kappa.score <- threeClassKappa(lead, scnd)
  kappa.score.backwards <- threeClassKappa(scnd, lead)
  stopifnot(kappa.score$global$kappa == kappa.score.backwards$global$kappa)
  ret.list <- list("procrustes.corr"=pseudo.procestes.corr, "factor.procrust.corr" = aligned.corr$factor.cors,"multi.class"=kappa.score, 
                   "rss"=rss.by.procrust, "kendall.corr"=aligned.corr.k$greedy.match$corr, "true.procrustes.corr"=cor.by.procrust)
  if(SWAP)
  {
    ret.list$order = "swapped";  ret.list$A = scnd;  ret.list[["B"]]=lead;  
    ret.list[["A.optimal"]] = aligned.corr$B;  ret.list[["B.optimal"]]=aligned.corr$A
    ret.list[["caret.report.A.first"]]=prd$caret.report.bkwd; ret.list[["caret.report.B.first"]]=prd$caret.report.fwd
    ret.list[["metrics.A.first"]] = prd$metrics.bwkd;  ret.list[["metrics.B.first"]]=prd$metrics.fwd 
    ret.list[["pr.A.first"]]=prd$pr.second; ret.list[["pr.B.first"]]=prd$pr.first
  }else
  {
    ret.list$order = "preserved";  ret.list$A = lead;  ret.list[["B"]]=scnd;
    ret.list[["A.ordered"]] =aligned.corr$A;  ret.list[["B.ordered"]]=aligned.corr$B
    ret.list[["caret.report.A.first"]]=prd$caret.report.fwd; ret.list[["caret.report.B.first"]]=prd$caret.report.bkwd
    ret.list[["metrics.A.first"]] =prd$metrics.fwd;  ret.list[["metrics.B.first"]]=prd$metrics.bwkd
    ret.list[["pr.A.first"]]=prd$pr.first; ret.list[["pr.B.first"]]=prd$pr.second
  }
  ret.list
  
}

#Get a binary PRcurve by looking at the columsn of the lead and second maftrix
#Binarized such that 0/1
getBinaryPRCurve <- function(lead, scnd)
{
  library(caret)
  lead.binarized <- lead
  scnd.binarized <- scnd
  lead.binarized[lead.binarized != 0] <- 1
  scnd.binarized[scnd.binarized != 0] <- 1
  predicted <- as.vector(scnd.binarized)
  actual <- as.vector(lead.binarized)
  #Using one first
  xtab <- table(predicted, actual)
  #If entries are all ones, than skip this
  if(all(predicted == 1) | all(actual == 1))
  {
    conf.forward <- NA;
    conf.backward <- NA
  }else
  {
    conf.forward <- caret::confusionMatrix(table(predicted, actual))
    conf.backward <- caret::confusionMatrix(table(actual, predicted))
  }
  
  #alternatively 
  library(Metrics)
  forward.metrics <- list("recall"=Metrics::recall(actual = actual, predicted = predicted),
                          "precision" = Metrics::precision(actual = actual, predicted = predicted),
                          "f1"=   Metrics::f1(actual, predicted))
  
  reverse.metrics <- list("recall"=Metrics::recall(actual = predicted, predicted = actual),
                          "precision" = Metrics::precision(actual = predicted, predicted = actual),
                          "f1"=   Metrics::f1(predicted, actual))
  
  #PR curve:
  scnd.scores <- scnd^2/max(scnd^2)
  #scores of all the trues.
  positives <- scnd.scores[lead.binarized != 0]
  negatives <- scnd.scores[lead.binarized == 0]
  curve.lead.first <- pr.curve(scores.class0 = positives, scores.class1 = negatives, curve = TRUE)
  
  
  lead.scores <- lead^2/max(lead)
  #scores of all the trues.
  positives <- lead.scores[scnd.binarized != 0]
  negatives <- lead.scores[scnd.binarized == 0]
  curve.scnd.first <- pr.curve(scores.class0 = positives, scores.class1 = negatives, curve = TRUE)
  
  
  return(list("caret.report.bkwd"=conf.backward,"caret.report.fwd"=conf.forward,
              "metrics.bwkd" = reverse.metrics, "metrics.fwd"=forward.metrics, "pr.second"=curve.scnd.first, "pr.first"=curve.lead.first))
  
}

splitThreeClasses <- function(v)
{
  lead.classes <- v
  lead.classes[lead.classes > 0] <- 1
  lead.classes[lead.classes < 0] <- 2
  lead.classes
}
#TODO: verify this is same regardless of order (it should be, since looking at similarity)
threeClassKappa <- function(lead, scnd, by_factor = TRUE)
{
  lead.classes <- splitThreeClasses(lead)
  scnd.classes <- splitThreeClasses(scnd)
  lead.vector <- as.vector(lead.classes) #(stacks columns)
  scnd.vector <- as.vector(scnd.classes)
  global.scores <- psych::cohen.kappa(data.frame(lead.vector, scnd.vector))
  message("Now by factor...")
  #Do it by column too
  #Be warned- this results in very sparse setups and so can yield error pars that are invalid.
  by.col = NA
  if(by_factor)
  {
    by.col <- lapply(1:ncol(lead.classes), function(i) psych::cohen.kappa(data.frame(lead.classes[,i], scnd.classes[,i])))
  }

  return(list("global"=global.scores, "by.factor"=by.col))
}
#Do we run across all the factorizations, or do we just specify the falls?
#Do across all, and get the stats in an output file.
#Or maybe, even better, write some functions to call from a markdown notebook where you actualy make the figures, eh?

#read in all the finngen and ukbb results into lists

loadGLEANERResults <- function(dir)
{
  finngen.list <- list()
  for(f in list.files(path = dir, pattern = "*_final_dat.RData"))
  {
    load(paste0(dir, f))
    name = gsub(x = f, replacement = "", pattern = "_final_dat.RData")
    finngen.list[[name]] <- ret
  }
  return(finngen.list)
}

loadGLEANERResultsLists <- function(dir_finngen, dir_ukbb)
{
  finngen.list <- loadGLEANERResults(dir_finngen)
  ukbb.list <- loadGLEANERResults(dir_ukbb)
  return(list("finngen" = finngen.list, "ukbb"=ukbb.list))
}

#Perform the model comparison on a given matrix pair from finngen.list and ukbb.list
#n : the arg name to look at
#BAD CODING PRACTICE- we add directly to the table in this function and then update it accordingly.
#compareFactPrecision(n, finngen.list, ukbb.list, comp.v,comp.u,global.df.V,global.df.U,full.df.V,full.df.U)
compareFactPrecision <- function(n, finngen.list, ukbb.list, comp.v,comp.u, global.df.V,global.df.U, overall.df.V, overall.df.U)
{
  #print(n)
  fg.v <- as.matrix(finngen.list[[n]]$V)
  uk.v <- as.matrix(ukbb.list[[n]]$V)
  comp.v[[n]] <- compareModelMatricesComprehensive(fg.v,uk.v)
  #statistics overall
  global.df.V <- cbind(global.df.V, 
                     c(comp.v[[n]]$metrics.A.first$recall,
                       comp.v[[n]]$metrics.A.first$precision,
                       comp.v[[n]]$procrustes.corr,
                       paste0(ncol(fg.v), "-",ncol(uk.v)),
                       comp.v[[n]]$pr.A.first$auc.integral,comp.v[[n]]$multi.class$global$kappa,
                       comp.v[[n]]$kendall.corr, comp.v[[n]]$true.procrustes.corr, comp.v[[n]]$rss))
  #Statistics by factor:
  overall.df.V <- rbind(overall.df.V,data.frame("study" = n, "Corr"=comp.v[[n]]$factor.procrust.corr, "Kappas"=sapply(comp.v[[n]]$multi.class$by.factor, function(x) x$kappa)))
  fg.u <- as.matrix(finngen.list[[n]]$U)
  uk.u <- as.matrix(ukbb.list[[n]]$U)
  comp.u[[n]] <- compareModelMatricesComprehensive(fg.u,uk.u,full.procrust = TRUE)
  global.df.U <- cbind(global.df.U, c((comp.u[[n]]$metrics.A.first$recall),
                                  comp.u[[n]]$metrics.A.first$precision, 
                                  comp.u[[n]]$procrustes.corr,paste0(ncol(fg.u), "-",ncol(uk.u)),
                                  comp.u[[n]]$pr.A.first$auc.integral,comp.u[[n]]$multi.class$global$kappa,
                                  comp.u[[n]]$kendall.cor,comp.u[[n]]$true.procrustes.corr, comp.u[[n]]$rss))
  
  
  overall.df.U <- rbind(overall.df.U,data.frame("study" = n, 
                                                "Corr"=comp.u[[n]]$factor.procrust.corr, 
                                                "Kappas"=sapply(comp.u[[n]]$multi.class$by.factor, function(x) x$kappa)))
  
  return(list(comp.v, global.df.V, comp.u, global.df.U,overall.df.V,overall.df.U))
}

##
#Pass in the lists of runs, and perform the comparison tests for each, and write the results out in an easy-to-read table.
#Note that this is not a fast function, the comparison step is kinda slow.
#Returns table with data for Recall, precision, correlation, K counts, and AUPRC
#makeTableOfComparisons(no.shrink.results$finngen, no.shrink.results$ukbb)
#Re-factor this code: give global scores as well as distributional scores across factors
#direct.gamma.comparisons <- makeTableOfComparisons(std_block.results$finngen, std_block.results$ukbb)
makeTableOfComparisons <- function(finngen.list, ukbb.list, by_order = FALSE)
{
  both <- intersect(names(ukbb.list), names(finngen.list))
  if(by_order)
  {
    both = 1:length(ukbb.list)
  }
  #drop BIC-sklearn_K-MAX_COVAR
  #both <- both[-which(both == "BIC-sklearn_K-MAX_COVAR")]
  comp.v <- list()
  comp.u <- list()
  comp.xhat <- list()
  #These get the global statistics, summarizing the entire factorization
  global.df.V <- NULL
  global.df.U <- NULL
  #These get the statistics by factor, allowing us to make boxplots
  full.df.V <- NULL
  full.df.U <- NULL
  for(n in both)
  {
    nfd <- compareFactPrecision(n, finngen.list, ukbb.list, comp.v,comp.u,global.df.V,global.df.U,full.df.V,full.df.U)
    comp.v<-nfd[[1]]; global.df.V<-nfd[[2]]; comp.u<-nfd[[3]]; global.df.U<-nfd[[4]];full.df.V <- nfd[[5]]; full.df.U<- nfd[[6]]
    #Comparison of X directly
    #Make sure the orders are the same.
    stopifnot(all(finngen.list[[n]]$trait.names == ukbb.list[[n]]$trait.names))
    stopifnot(all(finngen.list[[n]]$snp.ids == ukbb.list[[n]]$snp.ids))
    Xhat.fg <- finngen.list[[n]]$U %*% t(finngen.list[[n]]$V)
    Xhat.uk <- ukbb.list[[n]]$U %*% t(ukbb.list[[n]]$V)
    comp.xhat<- c(comp.xhat, norm(Xhat.fg - Xhat.uk, "F"))
    }
  
  rownames(global.df.V) <- c("recall", "precision", "correlation", "Fg_K-UK_K", "PRC", "GlobalKappa", "Kendall_correlation", "Procrustes_pearson", "Euclidian_dist")
  colnames(global.df.V)<- both
  
  rownames(global.df.U) <- c("recall", "precision", "correlation", "Fg_K-UK_K", "PRC","GlobalKappa","Kendall_correlation", "Procrustes_pearson", "Euclidian_dist")
  colnames(global.df.U)<- both
  
  pr.dat.v <- data.frame(t(global.df.V)) %>% tibble::rownames_to_column("rn") %>% 
    separate(rn, into = c("BIC", "Factors"), sep = "_K") %>% 
    mutate("BIC" = gsub(replacement = "", pattern = "BIC-", x = BIC),"Factors" = gsub(replacement = "", pattern = "-", x = Factors)) #%>% print()
  pr.dat.v$recall <- as.numeric(pr.dat.v$recall)
  pr.dat.v$precision <- as.numeric(pr.dat.v$precision)
  pr.dat.v$correlation <- as.numeric(pr.dat.v$correlation)
  pr.dat.v$PRC <- as.numeric(pr.dat.v$PRC)
  
   pr.dat.u <- data.frame(t(global.df.U)) %>% tibble::rownames_to_column("rn") %>% 
    separate(rn, into = c("BIC", "Factors"), sep = "_K") %>% 
    mutate("BIC" = gsub(replacement = "", pattern = "BIC-", x = BIC),"Factors" = gsub(replacement = "", pattern = "-", x = Factors)) #%>% print()
  pr.dat.u$recall <- as.numeric(pr.dat.u$recall)
  pr.dat.u$precision <- as.numeric(pr.dat.u$precision)
  pr.dat.u$correlation <- as.numeric(pr.dat.u$correlation)
  pr.dat.u$PRC <- as.numeric(pr.dat.u$PRC)
  
  #Directly compare the output X matrices. This may be a better metric
    #Add some more comprehensive data tables too
  return(list("V"=pr.dat.v, "U"=pr.dat.u, "X_hat.norms"=comp.xhat, "full.V" =full.df.V, "full.U" = full.df.U))
  
}

#' Compare all the covariance-adjusted outputs pairwise
#'
#' @param block.results a list containing two lists, one for finngen and one for ukbb, with the loaded RData in each entry
#'Note that here, the finngen ones are the ones changing each time, UKBB results stay the same.
#' @return list of lists containing the matching statistics for each pair gamma parameters
#' @export
#'
#' @examples
fullPairwiseComparisons <- function(block.results, len.run=11)
{
  reordered.lists.fg <- list()
  reordered.lists.uk <- list()
  stopifnot(length(block.results$finngen) == length(block.results$ukbb))
  adj.length = len.run-1
  for(i in 1:adj.length)
  {
    reordered.lists.fg[[i]] <- lapply(i:(adj.length+i), function(j) block.results$finngen[[(j%%len.run+1)]])
    names(reordered.lists.fg[[i]]) <- sapply(i:(adj.length+i), function(j) names(block.results$finngen)[(j%%len.run+1)])
    
    reordered.lists.uk[[i]] <- lapply(i:(adj.length+i), function(j) block.results$ukbb[[(j%%len.run+1)]])
    names(reordered.lists.uk[[i]]) <- sapply(i:(adj.length+i), function(j) names(block.results$ukbb)[(j%%len.run+1)])
  }
  #Need to also add the matching case to get the diagonal elements:
  reordered.lists.fg[[len.run]] <- block.results$finngen
  reordered.lists.uk[[len.run]] <- block.results$ukbb
  lapply(reordered.lists.fg, function(x) makeTableOfComparisons(block.results$ukbb, x,by_order = TRUE))
  #This "by_order" flag super important!
  
}


#' Function to take all the pairwise run data and put it into a table
#'
#' @param all.pairwise.settings a list of objects looking at the U and V comparison for pairiwise test
#' elements have U,V, Xhat norms, and more complete U and V data (distribution) [.e.g a list of objects outputted from makeTableOfComparisons]
#' @param block_results the results from both UKBB and finngen- the UKBB is used as the "unchanging" name reference list
#' @param changing.lists the lists of reordered gammas for a given sample, to match the "Finngen: reorderin, 
#' since that one had the rotating gammas in fullPairwiseComparisons above (different gamma references)
#' @param pattern the regualr expression of gamma terms to extract (group 1)
#'
#' @return
#' @export
#'
#' @examples
extractAndTabularizeRun <- function(all.pairwise.settings, block_results, changing.lists, pattern = "adjusted_([0-9\\.]+)-shrinkage")
{
  growing.V.crosscohort <- NULL
  growing.U.crosscohort <- NULL
  Xhat_norms <- NULL
  for(i in 1:length(all.pairwise.settings))
  {
    growing.V.crosscohort <- rbind(growing.V.crosscohort, all.pairwise.settings[[i]]$V %>% 
                                     mutate("UKBB.factors" = names(block_results$ukbb), "FG.factors"= names(changing.lists[[i]])) %>%
                                     mutate("Gamma.ukbb" = as.factor(seq(0,1,0.1))))
    growing.U.crosscohort <- rbind(growing.U.crosscohort, all.pairwise.settings[[i]]$U %>% 
                                     mutate("UKBB.factors" = names(block_results$ukbb), "FG.factors"= names(changing.lists[[i]]))%>%
                                     mutate("Gamma.ukbb" = as.factor(seq(0,1,0.1))))
    Xhat_norms <- rbind(Xhat_norms, data.frame("dist" = unlist(all.pairwise.settings[[i]]$X_hat.norms),
                                               "Gamma.fg"=as.factor(stringr::str_match(names(changing.lists[[i]]), pattern)[,2]),  "Gamma.ukbb"=as.factor(seq(0,1,0.1))))
  }
  growing.V.crosscohort$Gamma.fg <- as.factor(stringr::str_match(growing.V.crosscohort$FG.factors, pattern)[,2])
  growing.U.crosscohort$Gamma.fg <- as.factor(stringr::str_match(growing.U.crosscohort$FG.factors, pattern)[,2])
  #Make all the relevant columns numeric:
  growing.V.crosscohort <- growing.V.crosscohort %>%
    mutate_at(c("correlation","PRC","GlobalKappa", "Kendall_correlation","Procrustes_pearson","Euclidian_dist"), as.numeric)
  growing.U.crosscohort <- growing.U.crosscohort %>%
    mutate_at(c("correlation","PRC","GlobalKappa", "Kendall_correlation","Procrustes_pearson","Euclidian_dist"), as.numeric)
  return(list("V" = growing.V.crosscohort %>% select(-Factors), 
              "U"=growing.U.crosscohort %>% select(-Factors), "Xhat_norms"=Xhat_norms))
}


#' Function to take all the pairwise run data and put it into a table
#'
#' @param all.pairwise.settings a list of objects looking at the U and V comparison for pairiwise test
#' elements have U,V, Xhat norms, and more complete U and V data (distribution) [.e.g a list of objects outputted from makeTableOfComparisons]
#' @param block_results the results from both UKBB and finngen- the UKBB is used as the "unchanging" name reference list
#' @param changing.lists the lists of reordered gammas for a given sample, to match the "Finngen: reorderin, 
#' since that one had the rotating gammas in fullPairwiseComparisons above (different gamma references)
#' @param pattern the regualr expression of gamma terms to extract (group 1)
#'
#' @return
#' @export
#'
#' @examples
extractAndTabularizeRunGeneric <- function(all.pairwise.settings)
{
  growing.V.crosscohort <- NULL
  growing.U.crosscohort <- NULL
  Xhat_norms <- NULL
  for(i in 1:length(all.pairwise.settings))
  {
    growing.V.crosscohort <- rbind(growing.V.crosscohort, all.pairwise.settings[[i]]$V)
    growing.U.crosscohort <- rbind(growing.U.crosscohort, all.pairwise.settings[[i]]$U)
    Xhat_norms <- rbind(Xhat_norms, data.frame("dist" = unlist(all.pairwise.settings[[i]]$X_hat.norms)))
  }
  #Make all the relevant columns numeric:
  growing.V.crosscohort <- growing.V.crosscohort %>%
    mutate_at(c("correlation","PRC","GlobalKappa", "Kendall_correlation","Procrustes_pearson","Euclidian_dist"), as.numeric)
  growing.U.crosscohort <- growing.U.crosscohort %>%
    mutate_at(c("correlation","PRC","GlobalKappa", "Kendall_correlation","Procrustes_pearson","Euclidian_dist"), as.numeric)
  return(list("V" = growing.V.crosscohort, 
              "U"=growing.U.crosscohort, "Xhat_norms"=Xhat_norms))
}


