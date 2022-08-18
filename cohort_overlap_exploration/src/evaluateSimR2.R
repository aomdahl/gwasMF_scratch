#### Function to evaluate error between predicted L and F and true L and F
#### Taken from Yuan He's code with some modification.

#Changes:
  #Sig_hits not considered 
#lp <- pred.loadings
#fp <- pred.factors
#
#evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.factors)
if(FALSE)
{
  trueL=true.loadings
  trueF = true.factors
  lp= pred.loadings
  fp=pred.factors
}

evaluteFactorConstruction <- function(trueL, trueF, lp, fp, easy=TRUE, verbose = FALSE){
  if(is.null(lp) | length(lp) == 0)
  {
    message("empty table")
    lp <- matrix(0, nrow= nrow(trueL), ncol = nrow(trueL))
  }
  if(is.null(fp) | length(fp) == 0)
  {
    message("empty table")
    fp <- matrix(0, nrow= nrow(trueF), ncol = nrow(trueF))
  }
  
  #make sure all are matrices
  
  
  trueL <- as.matrix(trueL); trueF <- as.matrix(trueF); lp <- as.matrix(lp); fp <- as.matrix(fp)
  
  
  #fill in missing columns with 0s
  if(ncol(fp) < ncol(trueF)){
    message("missing factors, will set to 0...")
    dif = ncol(trueF) - ncol(fp)
    fp = cbind(fp, matrix(rep(0, nrow(fp) * dif), ncol = dif))
    lp = cbind(lp, matrix(rep(0, nrow(lp) * dif), ncol = dif))
  }
  #fill in NAs
  fp[is.na(fp)] = 0 
  lp[is.na(lp)] = 0
  
  K = ncol(trueF)
  suppressWarnings(library('combinat'))
  if(K > 9)
  {
    message("There are too many latent factors for comprehensive examination of options. Do a greedy search.")
    return()
  } else if(easy) 
  {
    ordering = permn(K)
    sign.ordering <- gtools::permutations(2, K, c(1,-1), repeats.allowed = TRUE)
    all.l <- getallR2(sign.ordering, ordering, lp, trueL)
    all.f <- getallR2(sign.ordering, ordering, fp, trueF)
    tot = all.l + all.f
    if(all(is.na(all.l)) & all(is.na(all.f)))
    {
      return(c(0,0))
    }
    return(c(all.l[which.max(tot)],all.f[which.max(tot)]))
    
  }
  else {
    message("Doing this the easy way...")
    #get all possible permutations in F and L
    #need to account for differences in possible sign assignments too....
    ordering = permn(K)
    f.perf <- tableCor(fp, trueF, ordering)
    l.perf <- tableCor(lp, trueL, ordering)
    #Look at the sum of correlations overall (taking absolute value)
    #Note that there may be an exceptional case, where they have opposite directions of correlation, and this is missed
    ord_sum = (f.perf$r2 + l.perf$r2) 
    ord = which.max(ord_sum)  #which combination is highest overall
    #Ensure the sign choice is the same - if one is a switched sign, the other must be too....
    if(!(all(f.perf$signs[[ord]] == l.perf$signs[[ord]])) & !(all(f.perf$signs[[ord]] == -1*l.perf$signs[[ord]]))) #also acceptable if they are just inverses of each other, that's an easy fix.
    {
      message("Sign assignments aren't consistent")
      message("Comparing an ordered approach...")
      #Here, fit L (done above), then use those orders on F
      start.pref <- list("l.first" =  tableCor(fp, trueF, ordering, sign.prior = l.perf$signs), 
                                "f.first" =  tableCor(lp, trueL, ordering, sign.prior = f.perf$signs))
      #Determine which of these yields the best 
      max.choices = (list("l.first" = max(l.perf$r2 + start.pref$l.first$r2), "f.first" = max(f.perf$r2 + start.pref$f.first$r2)))
      #1 is lfirst, 2 is f.first
      #Choose the maximizing one.
      max.fit = start.pref[[which.max(max.choices)]]
      #get the order index of the maximizing one
      ord.second <- which.max(max.fit$r2)
      if(names(which.max(max.choices)) == "l.first")
      {
        l_corr = l.perf$r2[[ord.second]]
        f_corr = start.pref$l.first$r2[[ord.second]]
      } else if (names(which.max(max.choices)) == "f.first"){
        f_corr = f.perf$r2[[ord.second]]
         l_corr = start.pref$f.first$r2[[ord.second]]
      }
      else
      {
        message("ERROR")
        l_corr = NA
        f_corr = NA
      }

    } else {
      l_corr = l.perf$r2[[ord]]
      f_corr = f.perf$r2[[ord]]
    }
    return(c(l_corr, f_corr))
  }
  
}

#note- there are 2 approaches here:
#1 we try each sign combo for each max combo. O(k^2*k!), very expensive for large K, but cheaper for small
#2 we first identify the best sign for each, then get the sums at the end (O(k^3 _ k!), still expensive but much cheaper

#function 
bestR2BySign <- function(pred.v, true) #takes in a full vector version.
{
  library(RcppAlgos)
  REF <- c(1,-1)
  K <- ncol(true)
  N <- nrow(true)
  true.v <- as.vector(true)
  #remove the last one, its all negatives
  sign.options <- RcppAlgos::permuteGeneral(REF, K, TRUE)
  r.choices <- c()
  for(i in 1:nrow(sign.options))
  {
    sign.expand <- unlist(lapply(sign.options[i,], function(x) rep(x, N)))
    r.choices <- c(r.choices, r2Signed(pred.v, sign.expand, true.v))
  }
  best <- which.max(r.choices)
  return(list("r2" = max(r.choices), "signs"= sign.options[best,]))
}

r2Signed <- function(pred.v, signs, true.v)
{
  return(cor(pred.v * signs, as.vector(true.v))^2)
}


#wrapper function to check all the possible orderings.
tableCor <- function(pred, true, ordering, sign.prior = NULL){
  f_cor = rep(0, length(ordering))
  f_signs = list()
  for(ord in seq(1, length(ordering))){
    curr.order <- as.vector(pred[,ordering[[ord]]]) #for the current ordering, which combination of signs yields the best fit?
    if(is.null(sign.prior))
    {
      best.sign <- bestR2BySign(curr.order, true)
    }else
    {
      best.sign <- list("r2"=r2Signed(curr.order, sign.prior[[ord]], true), "signs"=sign.prior[[ord]])
    }
    
    f_cor[ord] = best.sign$r2
    f_signs[[ord]]=best.sign$signs
  }
  return(list("r2" = f_cor, "signs"=f_signs))
}


getallR2 <- function(sign.ordering, ordering, pred, actual)
{
  nreps = nrow(actual)
  all.options.factors <- apply(sign.ordering, 1, function(s) {smat <- do.call("rbind", lapply(1:nreps, function(x) s)); lapply(ordering, function(x) pred[,x] * smat)})
  res.mat <- matrix(NA,nrow(sign.ordering),length(ordering))
  #For each of those, estimate the R2
  for(signopt in 1:nrow(sign.ordering))
  {
    res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) cor(as.vector(x), as.vector(actual))^2))
  }
  return(res.mat)
}



