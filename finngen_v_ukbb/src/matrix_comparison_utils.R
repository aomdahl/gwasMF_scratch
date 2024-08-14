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

#' Assessment of matrix similarity using several different types of metrics
#' Approaches include pearson and kendall correlation, classification error (Precision-recall curve),
#' categorization similarity (Cohen's Kappa metric), and Euclidian distance (frobenius norm)
#' All of these require columns to be mapped to each other, which can be done either by a greedy correlation search, or by a procrustes analysis
#' This function automatically does both
#' @param A 
#' @param B 
#' @param corr.type 
#' @param full.procrust 
#'
#' @return
#' @export
#'
#' @examples
compareModelMatricesComprehensive <- function(A,B, corr.type = "pearson")
{
  #First, align the matrices to have the same dimensions
  library(PRROC)
  prepped <- prepMatricesForAnalysis(A,B)
  lead <- prepped$lead;  scnd <- prepped$second; SWAP <- prepped$swap
  
  #Use my homeade correlation matching
  pseudo.procrust.results <- corrBasedScoring(lead, scnd, corr.method = "pearson")
  pseudo.ordered.matrices <- pseudo.procrust.results$corr_ordered_matrices
  
  #Use a published procrustes package too
  procrust.results <- procrustesBasedScoring(lead,scnd)
  procrust.ordered.matrices <- procrust.results$procrust_ordered_matrices

  #Need to verify and check orders are consistent here, before we do downstream analysis
  #checkUVAlignment
  #All we need to update here- lead and scnd, x_pearson, x_kendall, x_rss, 
  #  corr_based_stats <- with(comp.u[[n]], 

  #repairUVAlignment <- function(v.dat, u.dat, original.fg.v, original.fg.u,original.uk.v,original.uk.u, source="procrust")
  
  #assess with greedy correlation based order
  by_pseudo_procrust_order <- classificationScoring(pseudo.ordered.matrices$lead,pseudo.ordered.matrices$scnd, SWAP)
  by_true_procrust_order <- classificationScoring(procrust.ordered.matrices$lead, 
                                                  procrust.ordered.matrices$scnd, SWAP)
  names(by_true_procrust_order) <- gsub(names(by_true_procrust_order), pattern = "^", replacement = "procrust_")
  names(by_pseudo_procrust_order) <- gsub(names(by_pseudo_procrust_order), pattern = "^", replacement = "corr_based_")
  
  #Update results if all entries are 0:
  if((all(scnd == 0) & !all(lead == 0)) | (all(lead == 0) & !all(scnd == 0)))
  {
    message("An empty factorization appeared. Correlation scores being set to 0.")
    ret.list$procrustes.corr=0;ret.list$kendall.corr=0;ret.list$true.procrustes.corr=0
  }
  
  #Update the return data
  ret.list <- c(pseudo.procrust.results,procrust.results,by_pseudo_procrust_order,by_true_procrust_order )
  #How closely do our versions compare?
  ret.list$methods_mismatch <- compareProcrustesAndHomeade(pseudo.procrust.results$corr_based_data, procrust.results$procrustes_data,A,B,lead,scnd,DEBUG=FALSE)
  if(SWAP)
  {
    ret.list$order = "swapped";  ret.list$A = scnd;  ret.list[["B"]]=lead;  
  }else
  {
    ret.list$order = "preserved";  ret.list$A = lead;  ret.list[["B"]]=scnd;
  }
  ret.list
}


#' Put the larger matrix first and fill empty columns with 0s
#'
#' @param A 
#' @param B 
#'
#' @returns list with objects
#' "lead"- the larger matrix
#' "second"- the smaller matrix, filled with 0s
#' "swap"- if we had to reorder from input
#' @export
#'
#' @examples
prepMatricesForAnalysis <- function(A,B)
{
  SWAP = FALSE #refers to sizing order to get to the same size.
  lead <- as.matrix(A); scnd <- as.matrix(B)
  if(ncol(B) > ncol(A))
  {
    SWAP = TRUE
    lead <- as.matrix(B)
    scnd <- as.matrix(A)
  }
  ret <- fillWithZeros(lead, scnd)
  scnd <- ret$fp
  list("lead"=lead, "second"=scnd, "swap"=SWAP)
}

#' Perform scoring on 2 matrices based on classification
#'
#' @param lead the alleged reference matrix (larger of the 2), ORDERED
#' @param scnd the alleged predicted matrix, ORDERED
#' @param SWAP if initial A and B correspond to lead and scnd, this is FALSE
#'
#' @return
#' @export
#'
#' @examples
classificationScoring <- function(lead, scnd, SWAP)
{
  #do a simple BINARY precision variance assessment
  prd <-getBinaryPRCurve(lead, scnd)
  
  #Do a multi-class assessment, looking at +/0/-
  kappa.score <- threeClassKappa(lead, scnd)
  kappa.score.backwards <- threeClassKappa(scnd, lead)
  stopifnot(kappa.score$global$kappa == kappa.score.backwards$global$kappa)
  
  ret.list <- list("multi.class"=kappa.score)
  
  if(SWAP)
  {
    ret.list$order = "swapped";  
    ret.list[["A.optimal"]] = scnd;  ret.list[["B.optimal"]]=lead
    ret.list[["caret.report.A.first"]]=prd$caret.report.bkwd; ret.list[["caret.report.B.first"]]=prd$caret.report.fwd
    ret.list[["metrics.A.first"]] = prd$metrics.bwkd;  ret.list[["metrics.B.first"]]=prd$metrics.fwd 
    ret.list[["pr.A.first"]]=prd$pr.second; ret.list[["pr.B.first"]]=prd$pr.first
  }else
  {
    ret.list$order = "preserved";  
    ret.list[["A.ordered"]] =lead;  ret.list[["B.ordered"]]=scnd
    ret.list[["caret.report.A.first"]]=prd$caret.report.fwd; ret.list[["caret.report.B.first"]]=prd$caret.report.bkwd
    ret.list[["metrics.A.first"]] =prd$metrics.fwd;  ret.list[["metrics.B.first"]]=prd$metrics.bwkd
    ret.list[["pr.A.first"]]=prd$pr.first; ret.list[["pr.B.first"]]=prd$pr.second
  }
  
  ret.list
}

#' Title
#'
#' @param lead alleged reference matrix, unordered
#' @param scnd alleged predicted matrix, unordered
#'
#' @return
#' @export
#'
#' @examples
procrustesBasedScoring <- function(lead,scnd, DEBUG=FALSE)
{
  message("Doing full procrustes analysis")
  procestes.corr <- fullProcrustes(lead, scnd) #This is too slow.

    procrust.ordered.matrices <- list("lead"=lead, 
                                    "scnd"=matrixSignsProduct(scnd[,procestes.corr$mapping_reorder],procestes.corr$mapping_signs))
  cor.by.procrust <- stackAndAssessCorr(procestes.corr$A.normalized, procestes.corr$B.transformed) #correlation with the raw adjusted matrices
  rss.by.procrust <- procestes.corr$RSS

  corr.k.by.procrust <- stackAndAssessCorr(procestes.corr$A.normalized, procestes.corr$B.transformed, cor.type = "kendall")
  
  list("procrustes_data"=procestes.corr, "procrust_ordered_matrices" = procrust.ordered.matrices, 
       "procrustes_pearson"=cor.by.procrust, "procrustes_rss"=rss.by.procrust, "procrustes_kendall"=corr.k.by.procrust)

}

#' Correlation and RSS scoring based on greedy correlation matching
#'
#' @param lead alleged reference matrix, unordered
#' @param scnd alleged predicted matrix, unordered
#' @param corr.method - which correlation to use to order factors, default is pearson
#'
#' @return
#' @export
#'
#' @examples
corrBasedScoring <- function(lead, scnd, corr.method = "pearson")
{
  #G et the match and the greedy correlation corresponding to it
  #Pseudo-procrustes is based on correlation, rather than euclidian distance
  aligned.corr <- pseudoProcrustes(lead, scnd,corr.type = corr.method) #TODO:check this is reasonably optimal.
  pseudo.procestes.corr <- aligned.corr$greedy.match$corr
  corr_ordered_matrices = list("lead"=aligned.corr$A, "scnd"=aligned.corr$B)
  #alternative: kendall:
  #aligned.corr.k <- pseudoProcrustes(lead, scnd,corr.type = "kendall") #TODO:check this is reasonably optimal. #this is much harder to calculate for U, weird still slow
  #^this version reorders by kendall
  #Below- retain the pearson-based order
  aligned.corr.k <- stackAndAssessCorr(aligned.corr$A, aligned.corr$B, cor.type = "kendall") #This step is slow for some reason....
  
  RSS <- norm(aligned.corr$A - aligned.corr$B,  type = "F")
  
  list("corr_based_data"=aligned.corr, "corr_ordered_matrices" = corr_ordered_matrices, 
       "corr_based_pearson"= pseudo.procestes.corr, "corr_based_rss"=RSS, "corr_based_kendall"=aligned.corr.k)
}

compareProcrustesAndHomeade <- function(corr.alignment, procrustes.alignment,A,B,lead,scnd,DEBUG=TRUE)
{
  ret.dat <- FALSE
  #Check- how closely does my alignment match with the procrustes version
  non.matching <- sum(corr.alignment$greedy.match$order.pred != procrustes.alignment$mapping_reorder)
  k.diff <- abs(ncol(A) - ncol(B)) #possible mismatches due to missing columns
  if(non.matching > k.diff)
  {
    message("In this test, our greedy approach yields a different ordering than true procrustes.")
    message("Review closely all possible outputs")
    ret.dat <- TRUE
  }
  if(DEBUG) #Look more closely at the correlation differences, don't need anymore
  {
    #Reorder according to procrustes
    pred.one <- as.matrix(scnd[,procestes.corr$mapping_reorder]) %*% diag(procestes.corr$mapping_signs);      lead.one <- lead
    cor.reorder <- stackAndAssessCorr(lead.one, pred.one)
    #Try the different possible settings
    if(cor.reorder > pseudo.procestes.corr)
    {
      message("Procrustes orientation yields a better correlation than ours.")
    }else if(cor.reorder < pseudo.procestes.corr){
      message("Our version outperforms procrustes analysis! No change to make")
    }else
    {
      message("The correlations are the same;, no change to make.")
    }
  }
  ret.dat
}

#Get a binary PRcurve by looking at the columsn of the lead and second maftrix
#Binarized such that 0/1
getBinaryPRCurve <- function(lead, scnd)
{
  #BAD ashton- don't fix code like this. 
  #return(list("caret.report.bkwd"=NA,"caret.report.fwd"=NA,
  #            "metrics.bwkd" = NA, "metrics.fwd"=NA, "pr.second"=NA, "pr.first"=NA))
  library(caret)
  if(any(is.na(lead))|any(is.na(scnd)))
  {
    message("Found some NA, not good")
    print(lead)
    print(scnd)
  }
  lead.binarized <- lead
  scnd.binarized <- scnd
  lead.binarized[lead.binarized != 0] <- 1
  scnd.binarized[scnd.binarized != 0] <- 1
  predicted <- as.vector(scnd.binarized)
  actual <- as.vector(lead.binarized)
  #Using one first
  xtab <- table(predicted, actual)
  #print(xtab)
  #If entries are all ones, than skip this
  if(all(predicted == 1) | all(actual == 1) | nrow(xtab) < 2 | ncol(xtab) < 2 | all(predicted == 0) | all(actual == 0) )
  {
  #Special case to deal with
    if(all(predicted==1) & !all(actual ==1))
    {
      warning("Special case where all classified the same- introducing random entry to estimate precision")
    }
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
  if(!all(scnd == 0))
  {
    scnd.scores <- scnd^2/max(scnd^2)
  }
  else
  {
    empt.list <- list("type"="PR", "auc.integral"=NA,"auc.davis.goadrich"=NA)
    return(list("caret.report.bkwd"=conf.backward,"caret.report.fwd"=conf.forward,
                "metrics.bwkd" = reverse.metrics, "metrics.fwd"=forward.metrics, "pr.second"=empt.list, "pr.first"=empt.list))
  }
  
  
  #scores of all the trues.
  positives <- scnd.scores[lead.binarized != 0]
  negatives <- scnd.scores[lead.binarized == 0]
  curve.lead.first <- pr.curve(scores.class0 = positives, scores.class1 = negatives, curve = TRUE)
  
  
  lead.scores <- lead^2/max(lead^2)
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
threeClassKappa <- function(lead, scnd, by_factor = FALSE)
{
  lead.classes <- splitThreeClasses(lead)
  scnd.classes <- splitThreeClasses(scnd)
  lead.vector <- as.vector(lead.classes) #(stacks columns)
  scnd.vector <- as.vector(scnd.classes)
  global.scores <- suppressMessages(psych::cohen.kappa(data.frame(factor(lead.vector), factor(scnd.vector))))
  #message("Now by factor...")
  #Do it by column too
  #Be warned- this results in very sparse setups and so can yield error pars that are invalid.
  by.col = NA
  if(by_factor)
  {
    #BEWARE- suppressing warnings in this part.
    by.col <- lapply(1:ncol(lead.classes), function(i) suppressWarnings(psych::cohen.kappa(data.frame(lead.classes[,i], scnd.classes[,i]))))
  }

  return(list("global"=global.scores, "by.factor"=by.col))
}
#Do we run across all the factorizations, or do we just specify the falls?
#Do across all, and get the stats in an output file.
#Or maybe, even better, write some functions to call from a markdown notebook where you actualy make the figures, eh?

#read in all the finngen and ukbb results into lists

loadGLEANERResults <- function(dir, data_source)
{
  finngen.list <- list()
  for(f in list.files(path = dir, pattern = paste0("*", data_source)))
  {
    load(paste0(dir, f))
    name = gsub(x = f, replacement = "", pattern = data_source)
    if(data_source == "_final_dat.RData")
    {
      finngen.list[[name]] <- ret
    }else if(data_source == "_bic_dat.RData")
    {
      finngen.list[[name]] <- bic.dat
      finngen.list[[name]]$V <- finngen.list[[name]]$optimal.v
      finngen.list[[name]]$U <- finngen.list[[name]]$optimal.u
    }else
    {
      message("Invalid dataset. Returning")
      return(NA)
    }
    
  }
  return(finngen.list)
}

loadGLEANERResultsLists <- function(dir_finngen, dir_ukbb, dat_source="_final_dat.RData")
{
  finngen.list <- loadGLEANERResults(dir_finngen,dat_source)
  ukbb.list <- loadGLEANERResults(dir_ukbb,dat_source)
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
  cols.corr_based = c("cb_recall", "cb_precision", "cb_pearson", "cb_kendall", "cb_rss", "cb_auc", "cb_kappa")
  corr_based_stats <- with(comp.v[[n]], 
                           c(corr_based_metrics.A.first$recall,corr_based_metrics.A.first$precision,
                             corr_based_pearson,corr_based_kendall,corr_based_rss,
                             corr_based_pr.A.first$auc.integral,corr_based_multi.class$global$kappa))
  
  cols.procrustes = c("pr_recall", "pr_precision", "pr_pearson", "pr_kendall", "pr_rss", "pr_auc", "pr_kappa")
  procrust_based_stats <- with(comp.v[[n]], 
                           c(procrust_metrics.A.first$recall,procrust_metrics.A.first$precision,
                             procrustes_pearson,procrustes_kendall,procrustes_rss,
                             procrust_pr.A.first$auc.integral, procrust_multi.class$global$kappa)) 
  
  global.df.V <- cbind(global.df.V, c(paste0(ncol(fg.v), "-",ncol(uk.v)), corr_based_stats, procrust_based_stats))
  rownames(global.df.V) <- c("Fg_K-UK_K",cols.corr_based,cols.procrustes)
  
  #Statistics by factor (DEPRECATED)
  #overall.df.V <- rbind(overall.df.V,data.frame("study" = n, "Corr"=comp.v[[n]]$factor.procrust.corr, "Kappas"=sapply(comp.v[[n]]$multi.class$by.factor, function(x) x$kappa)))
  overall.df.V <- NA
  fg.u <- as.matrix(finngen.list[[n]]$U)
  uk.u <- as.matrix(ukbb.list[[n]]$U)
  comp.u[[n]] <- compareModelMatricesComprehensive(fg.u,uk.u)

  corr_based_stats <- with(comp.u[[n]], 
                           c(corr_based_metrics.A.first$recall,corr_based_metrics.A.first$precision,
                             corr_based_pearson,corr_based_kendall,corr_based_rss,
                             corr_based_pr.A.first$auc.integral,corr_based_multi.class$global$kappa))
  
  procrust_based_stats <- with(comp.u[[n]], 
                               c(procrust_metrics.A.first$recall,procrust_metrics.A.first$precision,
                                 procrustes_pearson,procrustes_kendall,procrustes_rss,
                                 procrust_pr.A.first$auc.integral, procrust_multi.class$global$kappa)) 
  
  global.df.U <- cbind(global.df.U, c(paste0(ncol(fg.u), "-",ncol(uk.u)), corr_based_stats, procrust_based_stats))
  rownames(global.df.U) <- c("Fg_K-UK_K",cols.corr_based,cols.procrustes)
    #Deprecated
  #overall.df.U <- rbind(overall.df.U,data.frame("study" = n, 
  #                                              "Corr"=comp.u[[n]]$factor.procrust.corr, 
  #                                              "Kappas"=sapply(comp.u[[n]]$multi.class$by.factor, function(x) x$kappa)))
  overall.df.U <- NA

  return(list(comp.v, global.df.V, comp.u, global.df.U,overall.df.V,overall.df.U,
              checkUVAlignment(comp.v[[n]],comp.u[[n]],fg.v,fg.u,uk.v,uk.u, source="procrust"), 
              checkUVAlignment(comp.v[[n]],comp.u[[n]],fg.v,fg.u,uk.v,uk.u, source="corr_based")))
}


checkUVAlignment <- function(v.dat, u.dat, fg.v, fg.u,uk.v,uk.u, source="procrust")
{
  if(source == "procrust")
  {
    mismatches <- sum(v.dat$procrustes_data$mapping_reorder != u.dat$procrustes_data$mapping_reorder)
  }
  if(source == "corr_based")
  {
    mismatches <- sum(v.dat$corr_based_data$greedy.match$order.pred != u.dat$corr_based_data$greedy.match$order.pred)
  }
  
  k.off <- abs(ncol(fg.v) - ncol(uk.v)) #possible mismatches with 0 columns
  if(mismatches > k.off)
  {
    warning("U and V weren't aligned in the same fashion, you have been warned.")
    if(source == "corr_based")
    {
      test.fix <- repairUVAlignment(v.dat, u.dat, fg.v, fg.u,uk.v,uk.u, source=source)
    }

    return(TRUE)
  }
  return(FALSE)
  
}


repairUVAlignment <- function(v.dat, u.dat, original.fg.v, original.fg.u,original.uk.v,original.uk.u, source="procrust")
{
  if(source == "procrust")
  {
    mismatches <- sum(v.dat$procrustes_data$mapping_reorder != u.dat$procrustes_data$mapping_reorder)
  }
  if(source == "corr_based")
  {
    mismatches <- sum(v.dat$corr_based_data$greedy.match$order.pred != u.dat$corr_based_data$greedy.match$order.pred)
  }
  k.off <- abs(ncol(original.fg.v) - ncol(original.uk.v))
  if(mismatches > k.off)
  {
    #an easy option- if the number of mismatched < 5, just use the other version....
    #update to have 0s in it, shoot.
    ret <- fillWithZeros(original.fg.v, original.uk.v)
    original.uk.v <- ret$fp
    ret <- fillWithZeros(original.fg.u, original.uk.u)
    original.uk.u <- ret$fp
    
    if(mismatches <= 5)
    {
      #first, identify the mismatched columns (based on index in updated matrix)
      mismatch.v <- which((v.dat$corr_based_data$greedy.match$order.pred != u.dat$corr_based_data$greedy.match$order.pred)) #index of mismatches in new space and in reference (fg)
      v.cols.off <- v.dat$corr_based_data$greedy.match$order.pred[mismatch.v] #get the columns #s in ORIGINAL space 
      mismatch.u <- which((u.dat$corr_based_data$greedy.match$order.pred !=v.dat$corr_based_data$greedy.match$order.pred)) #index of mismatches in new space and in reference (U)
      u.cols.off <- u.dat$corr_based_data$greedy.match$order.pred[mismatch.u] #original column #s in U
      stopifnot(all(u.cols.off %in% v.cols.off)); stopifnot(length(u.cols.off) == length(v.cols.off)) #The column numbers should be the same in both cases, just ordered differently
      #extract those same columns from U and V and make sure they are in the right order
      sorted.off <- sort(u.cols.off) #Sort them so they are order consistent; indices from the ORIGINAL,unsorted matrix that correspond to the sorted matrix 
      
      #Extract just the columns to re-score
      fg.mismatch.v <- original.fg.v[,mismatch.u]
      uk.mismatch.v <- original.uk.v[,sorted.off]
      fg.mismatch.u <- original.fg.u[,mismatch.u]
      uk.mismatch.u <- original.uk.u[,sorted.off]
      
      #Extract the columns that are already in a satisfactory order
      fg.match.v <- v.dat$corr_based_data$A[,-mismatch.v]
      uk.match.v <- v.dat$corr_based_data$B[,-mismatch.v]
      fg.match.u <- u.dat$corr_based_data$A[,-mismatch.u]
      uk.match.u <- u.dat$corr_based_data$B[,-mismatch.u]
      
      #Find the order that maximizes R2
      fixed.order <- evaluate_error(fg.mismatch.u, fg.mismatch.v, uk.mismatch.u, uk.mismatch.v)
      
      #Build the new matrices, with order
      uk.match.u <- cbind(uk.match.u,fixed.order$reordered_U)
      uk.match.v <- cbind(uk.match.v,fixed.order$reordered_V)
      
      fg.match.u <- cbind(fg.match.u,fg.mismatch.u)
      fg.match.v <- cbind(fg.match.v,fg.mismatch.v )
      
      #Get the new correlation:
      #update the order list to reflect the unification
      
      #Check that the order is now consistent.
      
      return(list("fg_v"=fg.match.v, "fg_u"=fg.match.u, "uk_v"=uk.match.v, "uk_u"=uk.match.u, "v_corr"=stackAndAssessCorr(fg.match.v,uk.match.v),"u_corr"= stackAndAssessCorr(fg.match.u,uk.match.u))) 
      
    }else
    {
      message("not implemented yet....")
      reorder.u <- getNewOrder(u.dat)
      reorder.v <- getNewOrder(v.dat)
      u.order.rss <- distanceByOrder(original.fg.u, original.fg.v, original.uk.u, original.uk.v, reorder.u)
      v.order.rss <- distanceByOrder(original.fg.u, original.fg.v, original.uk.u, original.uk.v, reorder.v)
      if(u.order < v.order)
      {
        v.dat <- updateProcrustOrder(v.dat,reorder.u )
      }else
      {
        u.dat <- updateProcrustOrder(u.dat,reorder.v)
      }
      return(list(v.dat, u.dat))
    }
  }
}

repairUVAlignmentProcrustes <- function(matrix.dat, matrix.list)
{
  reorder.u <- getNewOrder(u.dat)
  reorder.v <- getNewOrder(v.dat)
  u.order.rss <- distanceByOrder(original.fg.u, original.fg.v, original.uk.u, original.uk.v, reorder.u)
  v.order.rss <- distanceByOrder(original.fg.u, original.fg.v, original.uk.u, original.uk.v, reorder.v)
  if(u.order < v.order)
  {
    v.dat <- updateProcrustOrder(v.dat,reorder.u )
  }else
  {
    u.dat <- updateProcrustOrder(u.dat,reorder.v)
  }
  return(list(v.dat, u.dat))
}


 getNewOrder <- function(mat.dat)
 {
   #return the original matrix order from the object
   mat.dat$procrustes_data$mapping_reorder
 }
 
 distanceByOrder <- function(trueU, trueV, predU,predV,new.order)
 {
   rrmse(trueU,predU[,new.order]) + rrmse(trueV,predV[,new.order])
 }
 
 #All we need to update here- lead and scnd, x_pearson, x_kendall, x_rss, 
 #  library(PRROC)
 #prepped <- prepMatricesForAnalysis(A,B)
 #lead <- prepped$lead;  scnd <- prepped$second; SWAP <- prepped$swap
 
 #Use my homeade correlation matching
 #pseudo.procrust.results <- corrBasedScoring(lead, scnd, corr.method = "pearson")
 #pseudo.ordered.matrices <- pseudo.procrust.results$corr_ordered_matrices
 
 #Use a published procrustes package too
 #procrust.results <- procrustesBasedScoring(lead,scnd)
 #procrust.ordered.matrices <- procrust.results$procrust_ordered_matrices
 updateProcrustOrder <- function(obj, new.order)
 {
   update <- obj
   #update matrices
   update$procrust_A.ordered
   update$procrust_B.ordered
   update$procrust_ordered_matrices$lead
   #update order
   update$procrustes_data$mapping_reorder
   update$procrustes_data$mapping_signs
   #update corr (rho and pearson)
   update$procrustes_pearson
   update$procrustes_kendall
   #update RSS
   update$procrustes_rss <- "test"
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
  comp.mimsatched.procrust <- c()
  comp.mimsatched.corr_based <- c()
  #These get the global statistics, summarizing the entire factorization
  global.df.V <- NULL
  global.df.U <- NULL
  #These get the statistics by factor, allowing us to make boxplots
  full.df.V <- NULL
  full.df.U <- NULL
  missing.dat <- c()
  for(n in both)
  {
    print(n)
    if(length(ukbb.list[[n]]) == 0 | (length(finngen.list[[n]]) == 0))
    {
      missing.dat <- c(missing.dat, n)
      global.df.V <- cbind(global.df.V, rep(NA,nrow(global.df.V)))
      global.df.U <- cbind(global.df.U, rep(NA,nrow(global.df.U)))
      comp.xhat<- c(comp.xhat,NA)
      next
    }
    if(n == "BIC-std_K-MAX")
    {
      #message("P.")
    }
    nfd <- compareFactPrecision(n, finngen.list, ukbb.list, comp.v,comp.u,global.df.V,global.df.U,full.df.V,full.df.U)
    comp.v<-nfd[[1]]; global.df.V<-nfd[[2]]; comp.u<-nfd[[3]]; global.df.U<-nfd[[4]];full.df.V <- nfd[[5]]; full.df.U<- nfd[[6]]
    #Comparison of X directly
    #Make sure the orders are the same.
    stopifnot(all(finngen.list[[n]]$trait.names == ukbb.list[[n]]$trait.names))
    stopifnot(all(finngen.list[[n]]$snp.ids == ukbb.list[[n]]$snp.ids))
    Xhat.fg <- finngen.list[[n]]$U %*% t(finngen.list[[n]]$V)
    Xhat.uk <- ukbb.list[[n]]$U %*% t(ukbb.list[[n]]$V)
    if(any(dim(Xhat.uk) != dim(Xhat.fg)))
    {
      message("Unsolved problem here... why?")
    }
    #updating to be the root mean squared error:
    #comp.xhat<- c(comp.xhat, norm(Xhat.fg - Xhat.uk, "F"))
    comp.xhat<- c(comp.xhat, sqrt((norm(Xhat.fg - Xhat.uk, "F")^2)/ (ncol(Xhat.fg) * nrow(Xhat.fg))))
    comp.mimsatched.procrust <-  c(comp.mimsatched.procrust, nfd[[7]])
    comp.mimsatched.corr_based <-  c(comp.mimsatched.corr_based, nfd[[8]])
    }
  
  missing_i <- which(both %in% missing.dat)
  #rownames(global.df.V) <- c("recall", "precision", "correlation", "Fg_K-UK_K", "PRC", "GlobalKappa", "Kendall_correlation", "Procrustes_pearson", "Euclidian_dist")
  colnames(global.df.V)<- both#[-missing_i]
  
  #rownames(global.df.U) <- c("recall", "precision", "correlation", "Fg_K-UK_K", "PRC","GlobalKappa","Kendall_correlation", "Procrustes_pearson", "Euclidian_dist")
  colnames(global.df.U)<- both#[-missing_i]
  
  pr.dat.v <- data.frame(t(global.df.V)) %>% tibble::rownames_to_column("rn") %>% 
    separate(rn, into = c("BIC", "Factors"), sep = "_K") %>% 
    mutate("BIC" = gsub(replacement = "", pattern = "BIC-", x = BIC),"Factors" = gsub(replacement = "", pattern = "-", x = Factors)) %>%
    mutate_at(colnames(.)[-c(1,2,3)], as.numeric)
  #pr.dat.v$recall <- as.numeric(pr.dat.v$recall)
  #pr.dat.v$precision <- as.numeric(pr.dat.v$precision)
  #pr.dat.v$correlation <- as.numeric(pr.dat.v$correlation)
  #pr.dat.v$PRC <- as.numeric(pr.dat.v$PRC)
  
   pr.dat.u <- data.frame(t(global.df.U)) %>% tibble::rownames_to_column("rn") %>% 
    separate(rn, into = c("BIC", "Factors"), sep = "_K") %>% 
    mutate("BIC" = gsub(replacement = "", pattern = "BIC-", x = BIC),"Factors" = gsub(replacement = "", pattern = "-", x = Factors))  %>%
     mutate_at(colnames(.)[-c(1,2,3)], as.numeric)
   #%>% print()
  #pr.dat.u$recall <- as.numeric(pr.dat.u$recall)
  #pr.dat.u$precision <- as.numeric(pr.dat.u$precision)
  #pr.dat.u$correlation <- as.numeric(pr.dat.u$correlation)
  #pr.dat.u$PRC <- as.numeric(pr.dat.u$PRC)
  
  #Directly compare the output X matrices. This may be a better metric
    #Add some more comprehensive data tables too
  return(list("V"=pr.dat.v, "U"=pr.dat.u, "X_hat.norms"=comp.xhat,"U-V-aligned_procrust"=comp.mimsatched.procrust,"U-V-aligned_corr_based"= comp.mimsatched.corr_based,
              "full.V" =full.df.V, "full.U" = full.df.U))
  
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


