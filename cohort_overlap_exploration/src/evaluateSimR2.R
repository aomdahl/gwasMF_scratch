#### Function to evaluate error between predicted L and F and true L and F
#### Taken from Yuan He's code with some modification.
unitNorm <- function(x)
{
  x/norm(x, "2")
}
unitNorms <- function(M)
{
  ret <- apply(M, 2, function(x) unitNorm(x))
  ret[is.na(ret)] <- 0
  ret
}

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
#as.matrix(true.v), as.matrix(FactorM)
#evaluateSingleMatrixConstruction(as.matrix(true.v), as.matrix(ret$V))
evaluateSingleMatrixConstruction <- function(truem, predm)
{
  truem <- as.matrix(truem); predm <- as.matrix(predm);

  #fill in missing columns with 0s
  if(ncol(predm) < ncol(truem)){
    #message("missing factors, will set to 0...")
    dif = ncol(truem) - ncol(predm)
    predm = cbind(predm, matrix(rep(0, nrow(predm) * dif), ncol = dif))
  }
  #What if the new one learns
  if(ncol(truem) < ncol(predm)){
    #message("missing factors, will set to 0...")
    dif = ncol(predm) - ncol(truem)
    truem = cbind(truem, matrix(rep(0, nrow(predm) * dif), ncol = dif))
  }


  #fill in NAs
  predm[is.na(predm)] = 0

  K = ncol(truem) #previously was truem, had to change because
  suppressWarnings(library('combinat'))
  if(K > 9)
  {
    message("There are too many latent factors for comprehensive examination of options. Do a greedy search.")
    return()
  } else
  {
    ordering = permn(K)
    sign.ordering <- gtools::permutations(2, K, c(1,-1), repeats.allowed = TRUE)
    all.l <- getallR2(sign.ordering, ordering, predm, truem)
    return(all.l[which.max(all.l)])
  }
}

#updated.mat <-fillWithZeros(trueF, fp, lp=lp)
fillWithZeros <- function(trueF, fp, lp = NULL)
{
  trueF <- as.matrix(trueF)
  fp <- as.matrix(fp)
  if(ncol(fp) < ncol(trueF)){
    #message("missing factors, will set to 0...")
    dif = ncol(trueF) - ncol(fp)
    fp = cbind(fp, matrix(rep(0, nrow(fp) * dif), ncol = dif))
    if(!is.null(lp))
    {
      lp = cbind(lp, matrix(rep(0, nrow(lp) * dif), ncol = dif))
    }
    
  }
  return(list("fp" = fp, "lp" = lp))
}

#Ordering may be specified as a list of lists, with possible orders
#sel.u, sel.v, res$U, res$V
#evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.factors,unit.scale = FALSE)
evaluteFactorConstruction <- function(trueL, trueF, lp, fp, easy=TRUE, verbose = FALSE, ordering = NULL,corr.type = "pearson",...){
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
  updated.mat <-fillWithZeros(trueF, fp, lp=lp)
  fp <- updated.mat$fp; lp <- updated.mat$lp

  #fill in NAs
  fp[is.na(fp)] = 0
  lp[is.na(lp)] = 0

  K = ncol(trueF)
  suppressWarnings(library('combinat'))
  if(K > 9 & is.null(ordering))
  {
    message("There are too many latent factors for comprehensive examination of options. Reverting to a greedy search, in which U and V are treated separately.")
    return(c(greedyMaxCorr(trueL, lp,cor.type = corr.type), greedyMaxCorr(trueF, fp,cor.type = corr.type)))
  } else if(easy)
  {
    if(is.null(ordering))
    {
      ordering = permn(K)
    }

    sign.ordering <- gtools::permutations(2, K, c(1,-1), repeats.allowed = TRUE)
    all.l <- getallR2(sign.ordering, ordering, lp, trueL,...)
    all.f <- getallR2(sign.ordering, ordering, fp, trueF,...)
    tot = all.l + all.f
    if(all(is.na(all.l)) | all(is.na(all.f)))
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

#getallR2(sign.ordering, ordering, predm, truem)
getallR2 <- function(sign.ordering, ordering, pred, actual, unit.scale = FALSE)
{
  nreps = nrow(actual)
  all.options.factors <- apply(sign.ordering, 1, function(s) 
    {smat <- do.call("rbind", lapply(1:nreps, function(x) s)); lapply(ordering, function(x) pred[,x] * smat)})
  res.mat <- matrix(NA,nrow(sign.ordering),length(ordering))
  #For each of those, estimate the R2
  if(unit.scale){actual <- unitNorms(actual) }
  for(signopt in 1:nrow(sign.ordering))
  {
    if(unit.scale)
    {
      
      res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) cor(as.vector(unitNorms(x)), as.vector(actual))^2))
       
    }else
    {
      res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) cor(as.vector(x), as.vector(actual))^2))
    }
    
  }
  return(res.mat)
}


#####5/31
###A greedy implementation to evaluate factor construction, one matrix at a time.
#The main recursive call- requires
#Correlation matrix between remaining factors
#R2 matrix of remaining factors
#signs- vector of the signs assigned to each column
#pairs-a list where each entry is a tuple of (reference vect, pred vect) that match to each other
#returns signs and pairs at the end.
#recurseGreedyMap(r.mat,r2, c(dir), list(dt))
recurseGreedyMap <- function(cor, r2, signs, pairs)
{
  dt <- which(r2 == max(r2), arr.ind = TRUE)[1,] #just go with the first one...
  dir <- sign(cor[dt[1], dt[2]])
  if(all(r2 == 0) | ncol(r2) == 1)
  {
    return(list("signs"=signs, "pairs" = pairs))
  }
  pairs[[length(pairs) + 1]] <- dt
  
  cor[dt[1],] <- 0; cor[,dt[2]] <- 0
  r2[dt[1],] <- 0; r2[,dt[2]] <- 0
  recurseGreedyMap(cor,r2, c(signs, dir), pairs)
}

#' Correlation between matrix columns optimized for my application
#'
#' @param x first matrix 
#' @param y second matrix
#' @param cor.type desired type (pearson, kendall currentlty accepted)
#'
#' @return matrix of correlation estimates (no p-values)
#' @export
customCorr <- function(x,y, cor.type = "pearson")
{
  if(cor.type == "kendall")
  {
    r.mat <- t(apply(x, 2, function(a) apply(y, 2, function(b) pcaPP::cor.fk(a,b))))
  }else
  {
    r.mat <- cor(x,y, method=cor.type);
  }
  r.mat
}
#Only does it for one vector, the one you specify first
#greedyMaxCorr(true.v, pred.f)
greedyMaxCorr <- function(true.v, fp, verbose = FALSE, cor.type = "pearson")
{
  ret <- fillWithZeros(true.v, fp)
  r.mat <- customCorr(true.v, ret$fp, cor.type=cor.type);

  #fast way to do this with cor.fk
  
  r.mat[is.na(r.mat)] <- 0
  #r2 <-  cor(true.v, ret$fp,method=cor.type)^2; r2[is.na(r2)] <- 0
  r2 <-  r.mat^2; r2[is.na(r2)] <- 0
  dt <- which(r2 == max(r2), arr.ind = TRUE)[1,]
  dir <- sign(r.mat[dt[1], dt[2]])
  
  #case- there are multiple maximums. Just pick one to start at.....
  
  #Zero-out the R terms that aren't helping us.....
  r.mat[dt[1],] <- 0; r.mat[,dt[2]] <- 0
  r2[dt[1],] <- 0; r2[,dt[2]] <- 0
  matched.cols <- recurseGreedyMap(r.mat,r2, c(dir), list(dt))

  order.list <- 1:ncol(true.v)
  true.order <- sapply(matched.cols$pairs, function(x) x[1]) %>%  c(., order.list[!(1:ncol(true.v) %in% .)])
  pred.order <- sapply(matched.cols$pairs, function(x) x[2])%>%  c(., order.list[!(1:ncol(true.v) %in% .)])
  pred.signs <-  c(matched.cols$signs, rep(1,(ncol(true.v) - length(matched.cols$signs))))
  #scale
  
  #pred.vector <- as.vector(as.matrix(ret$fp[,..pred.order] %>% unitNorms(.)) %*% diag(pred.signs))
  #true.vector <- as.vector(as.matrix(true.v[,..true.order] %>% unitNorms(.)))
  #cor(pred.vector,true.vector)
  pred.one <- as.matrix(ret$fp[,pred.order] %>% unitNorms(.)) %*% diag(pred.signs)
  pred.two <- (as.matrix(true.v[,true.order]) %>% unitNorms(.))
  c_ <- stackAndAssessCorr(pred.one,pred.two)
  if(!verbose)
  {
    return(c_)
  }else
  {
    return(list("corr"=c_, "order.true"=true.order, "order.pred"=pred.order, "signs"=pred.signs))
  }
  
}

testGreedyMaxCorr <- function()
{
  #Manual example for checking performacne
  A <- matrix(rnorm(50),nrow=10)
  B <- matrix(0,nrow=10, ncol = 5)
  for(i in 1:5)
  {
    print(i%%5+1)
    B[,(i%%5+1)] <- A[,i]*(i/3) + rnorm(10,sd=0.1) + rnorm(10,sd=0.03)
    
  }
  cor(A,B)
  #Correct order: 
  test <- greedyMaxCorr(A, B, verbose = TRUE, cor.type = "pearson")
  test.tab <- data.frame("A_order"= test$order.true, "B_order" =test$order.pred) %>% arrange(`A_order`)
  stopifnot(all(test.tab$B_order==c(2,3,4,5,1)))
  
  A.c <- A
  A.c[,1] <- A[,2]
  A.c[,2] <- A[,5]
  A.c[,3] <- A[,4]
  A.c[,4] <- A[,3]
  A.c[,5] <- A[,1]
  test2 <- greedyMaxCorr(A.c, B, verbose = TRUE, cor.type = "pearson")
  test2.tab <- data.frame("A_order"= test2$order.true, "B_order" =test2$order.pred) %>% arrange(`A_order`)
  stopifnot(all(test2.tab$B_order==c(3,1,5,4,2)))
  #Second test- add some zeroes in there.
  B[,3] <- rep(0,10)
  #Correct order: 
  test <- greedyMaxCorr(A, B, verbose = TRUE, cor.type = "pearson")
  test.tab <- data.frame("A_order"= test$order.true, "B_order" =test$order.pred) %>% arrange(`A_order`)
  stopifnot(all(test.tab$B_order==c(2,3,4,5,1)))
  
  #Shuffle A a bit up

  
}

pseudoProcrustes <- function(A,B, corr.type="pearson")
{
  ret <- fillWithZeros(A, B)
  #TODO- why do we change the order of both matrices? Dpesn't make sense to me....
  cor.dat <- greedyMaxCorr(A, ret$fp, verbose = TRUE, cor.type = corr.type)
  new.A <- as.matrix(A[,cor.dat$order.true]); new.B <- as.matrix(ret$fp[,cor.dat$order.pred] %*% diag(cor.dat$signs))
  list("A" = new.A, "B" = new.B, "greedy.match"=cor.dat, "factor.cors" = diag(customCorr(new.A, new.B, cor.type = corr.type)))
}

stackAndAssessCorr <- function(A,B, cor.type = "pearson")
{
  pred.vector <- as.vector(A)
  true.vector <- as.vector(B)
  customCorr(pred.vector,true.vector, cor.type = cor.type)
}

#Taken from https://rpubs.com/mengxu/procrustes_analysis
#procestes.corr <- fullProcrustes(lead, scnd)
procrustes <- function(A, B){
  
  if(ncol(A) ==1 & ncol(B)==1)
  {
    A.normalized <- scale(A)
    B.normalized <- scale(B)
    RSS <- norm(A.normalized - B.normalized,  type = "F")
    return(list(A.normalized = A.normalized, B.normalized = B.normalized, rotation.mtx = 1, B.transformed = B.normalized, RSS = RSS))
  }
  
  # center and normalize A 
  A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A.centered, type = "F") / (ncol(A) * nrow(A))
  A.normalized <- A.centered / A.size
  
  # center and normalize B
  B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B.centered, type = "F") / (ncol(B) * nrow(B))
  B.normalized <- B.centered / B.size

  # Rotation matrix T 
  #if the matrix is big, just get the first few SVs
  m = B.normalized %*% t(A.normalized)
  if(nrow(m) > 1000 & ncol(m) > 1000)
  {
    message("Large matrix, need to use fast SVD. Please wait...")
    #svd.results <- corpcor::fast.svd(m)
    #Executive decision to just do 20 for speed
    svd.results <- RSpectra::svds(m, 20)
  }else{
    svd.results <- svd(m)
  }
  
  U <- svd.results$u
  V <- svd.results$v
  tr <- V %*% t(U)
  
  # B transformed
  B.transformed <- tr %*% B.normalized
  
  # Error after superimposition
  RSS <- norm(A.normalized - B.transformed,  type = "F")
  
  # Return
  return(list(A.normalized = A.normalized, B.normalized = B.normalized, rotation.mtx = T, B.transformed = B.transformed, RSS = RSS))
}

# fullProcrustes(lead, scnd)
fullProcrustes <- function(A,B)
{
  ret <- fillWithZeros(A, B)
  procrustes(as.matrix(A), as.matrix(ret$fp))
  
}

stackAndAssessNaive <- function(A,B)
{
  ret <- fillWithZeros(A, B)
  stackAndAssessCorr(A, ret$fp)
}


