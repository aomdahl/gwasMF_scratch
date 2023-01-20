#Temp script with functions for running the BIC-based fitting approach.


reportSimilarities <- function(f, s, t)
{
  print("Alphas:")
  print(c(f$alpha, s$alpha, t$alpha))
  print("Lambdas:")
  print(c(f$lambda, s$lambda, t$lambda))
  
  print("Ks:")
  print(c(f$K, s$K, t$K))
}


#Estimator of the L1 degrees of freedom
MatrixDFV <- function(mat_in, fixed_first = FALSE){
  print("RED ALERT- USING THE BAD V")
  if(fixed_first)
  {
    sum(mat_in[,-1] == 0)
  }else
  {
    sum(mat_in == 0)
  }
 
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

#CalcMatrixBIC(X,W,U,V,df=df,fixed_first = fixed_first,...)
CalcMatrixBIC <- function(X,W,U,V, ev="std", weighted = FALSE, df = NULL, fixed_first = FALSE)
{
  #Calculated based on my understanding of what is given by equations 12 and 7 in Lee et al 2010
  #We assume U (First term) is known/fixed, V is our predictor we evaluate here
  #11/17 UPDATE: df is calculated for U and V!
  
  n=nrow(X)
  d = ncol(X)
  if(fixed_first)
  {
	  message("Removing first factor because its fixed")
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
    message("not weighting here..")
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
  { #message("in correct place...")
    #print(ev)
    errvar = ev 
    }

  if(is.null(df))
  {
	  message("something went wrong,..")
    df <- MatrixDFU(V)
  }
  if(FALSE)
  {
    ret <- norm(X*W - (U %*% t(V))*W, type = "F")^2/(n*d * errvar) + (log(n*d*k)/(n*d))*df
  }else
  {
   
    if(df == 0)
    {
      message("Stop here plz.")
      message("TODO- get the program to stop or reset or something.")
    }
    #message("First term ", norm(X*W - (U %*% t(V))*W, type = "F")^2/(n*d * errvar))
    #message("second term ", (log(n*d)/(n*d))*df)
  ret <- norm(X*W - (U %*% t(V))*W, type = "F")^2/(n*d * errvar) + (log(n*d)/(n*d))*df
  #message("Total ", ret)
 #message(". ")
  }
  
  return(ret)
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

#WeightedAltVar(t(X),t(W),initV)
#Determine the variance using OLS. Basically, this is the "best" residual variance we can get.

WeightColumnHelper <- function(w,x,U)
{
  return(list("wx"=w*x, "wu" = w*U))
}

RegressOutMatrixComponent <- function(X,W,mc)
{
  ret <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
  for(col in 1:ncol(X))
  {
    w <- unlist(W[,col])
    wx <- (w*X[,col])
    wu <- w * mc
    fit <- lm(wx~wu + 0)
    fitted <- mc*coef(fit)
    ret[,col] <- X[,col] - fitted
  }
  return(ret)
}

WeightedAltVar <- function(X,W,U)
{
  n <- nrow(X) * ncol(X)
  p <- ncol(U) * ncol(X)
  #after discuussion with eric, just trying mle
  #p <- 0
  #message("N:", n)
  #message("P:", p)
  resids <- c()
  for(col in 1:ncol(X))
  {
    w <- unlist(W[,col])
    wx <- (w*X[,col])
    wu <- w * U
    fit <- lm(wx~wu + 0)
    resids <- c(resids, resid(fit))
  }
  r <- var(resids)
  print("Default var is ")
  print(r)
  # residual variance is 1/(n-p) e^Te
  r = 0
  r = (1/(n-p))*sum(resids * resids)
  print("my calculated r is:")
  print(r)
  #This is equivilant to all ther terms being put into one matrix
  if((sum(resids * resids) == 0))
  {
    message("No residual variance with current fit.")
    print(r)
    message("Be wary here- not clear what to do. Going to just approximate to small number")
    message("This is strong evidence of overfitting- we recommend dropping current factors.")
    message("Current U dims:")
    print(dim(U))
    #this would be easy enough to do, if we had everything in "PVE" order.
    return(1e-10)
  }
  return(r)
}

#Helper calls
CalcVBIC <- function(X,W,U,V,fixed_first=FALSE,...)
{
  #Looking at the effect of a small change..
  #
  #CalcMatrixBIC(X,W,U,V,df=MatrixDFV(V),...) #here, we look at sparsity, not number of non-zero coefs.
  #df = MatrixDFV(V, fixed_first = fixed_first)
  df=MatrixDFU(V,fixed_first=fixed_first)
  #print(paste0("DF is ", df))
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
   #For the BIC calculation, we don't  want to include dense factor 1....
  #if(option$fixed_ubiq)
  #{
  #      message("please be sure to remove F1")
  #}

  
  
  if(weighted) {
    #message("Weighted way...")
    if(option$fixed_ubiq)
    {
      #If we are down to just 1 column, don't calculate BICs, there is no point anymore.
      if(ncol(initU) == 1)
      {
        #If we have fixed_ubiq and down to 1 column, it doesn't matter 
        return(list("fits" = f.fits, "BIC"=rep(0,length(lambdas))))
      }
      Xr <- RegressOutMatrixComponent(X,W,initU[,1])
      av <- WeightedAltVar(Xr,W,as.matrix(initU[,-1]))
    }else
    {
      av <- WeightedAltVar(X,W,initU)
    }
    
    bics <- unlist(lapply(f.fits, function(x) CalcVBIC(X,W,initU,x$V, ev=av, weighted = TRUE, fixed_first = option$fixed_ubiq)))
  }else{
    av <- AltVar(X,initU)
    bics <- unlist(lapply(f.fits, function(x) CalcVBIC(X,W,initU, x$V, ev=av, fixed_first = option$fixed_ubiq)))
  }
  
  if(Inf %in% bics)
  {
    message("Error here...")
  }
  return(list("fits" = f.fits, "BIC"=bics))
}

DropEmptyColumns <- function(matin)
{
  matin <- as.matrix(matin)
  c <- colSums(matin != 0)
  if(any(c == 0))
  {
    #print("dropping")
    drops <- which(c == 0)
    if(1 %in% drops){
      message("BEWARE- 1st column dropping???")
    }
    return(as.matrix(matin[,-drops]))
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
    #OLS remove first X
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
#@param curr.mode= the mode of the sparsity space based on the current V and U settings
#@return a list of new sparsity points to try out.


ProposeNewSparsityParams <- function(bic.list,sparsity.params, curr.dist, n.points = 7, no.score = FALSE, one.SD.rule = FALSE)
{
  curr.mode = DensityMode(curr.dist)
  global.min <- min(curr.dist)
  if(length(bic.list) == 1)
  {
    message("No list to choose from- have zeroed all out..?")
    #Go from cuyrrent value to the mode, give a spread
    return(sort(10^seq(log10(sparsity.params),log10(curr.mode),length.out=n.points)))
  }
  if(no.score)
  {
    #then bic.list is the optimal one; generate fake scores
    fake.scores <- rep(100,length(sparsity.params))
    fake.scores[which(sparsity.params == bic.list)] <- -1
    bic.list <- fake.scores
  }
  optimal.index <- selectOptimalScoreIndex(bic.list, sparsity.params, one.SD.rule)
  #cases with redundancy are complex.
  optimal.sparsity.param <- sparsity.params[optimal.index]
  sorted.sparsity.params <- sort(sparsity.params, index.return = TRUE)
  ordered.list <- sorted.sparsity.params$x
  sorted.index <- which(sorted.sparsity.params$ix == optimal.index)
  #what is the index int eh sorted list of my optimal sparsty parameter? 
  #New paradigm: always look above and below,
  #If its the smallest paramter tested
  if(min(ordered.list) == optimal.sparsity.param)
  {
    message('best case is the minimum..')
    above <- ordered.list[sorted.index + 1]
    #below <- optimal.sparsity.param - (above - optimal.sparsity.param)
    #simplify this: we are stil searching, so look orders of magnitude
    #below <- 1e-10
    below <- global.min
    if(below > above)
    {
      #message("Global min param of distribution is greater than current one")
      #print(above)
      #print(below)
      below <- optimal.sparsity.param/2
    }
  } else if(max(ordered.list) == optimal.sparsity.param) #its the largest parameter tested
  {
    message('best case is the maximum')
    below <- ordered.list[sorted.index - 1]
    above <- curr.mode
    #new.list <- 10^seq(log10(below),log10(above),length.out=n.points)
  }else {
    #Its bounded- our estimates should be between the one immediately above and below
    above <- ordered.list[sorted.index + 1]
    below <- ordered.list[sorted.index - 1]
    #new.list <- seq(below,above,length.out = n.points)
  }
  if(length(above) > 1 | length(below) < 1)
  {
    message("WHAT is going on...")
    print(above)
    print(below)
    print(bic.list)
    print(sparsity.params)
    readline()
    quit()
  }
  if(is.na(above) | is.na(below))
  {
    print("here")
  }
  if(above == below)
  {
    message("Converged on single solution")
    return(below)
  }
  new.list <- 10^seq(log10(below),log10(above),length.out=n.points)
  
  #Ensure none of them are less than 0
  
  
  if(any(new.list <= 0))
  {
    rep.list <- new.list[which(new.list > 0)]
    if(any(new.list == 0))
    {
      rep.list <- c(min(rep.list)/2,rep.list)
    }
    new.list <- rep.list
  }

  unique(sort(c(optimal.sparsity.param, new.list)))
}

quickSort <- function(tab, col = 1)
{
  tab[order(tab[,..col], decreasing = TRUE),]
}
    oneSDRule <- function(bics, params)
    {
      if(Inf %in% params)
      {
        message("issue here....")
      }
      if(length(unique(bics)) == 1)
      {
        #all the parameters are the same
        message("BIC score for all parameters are same. Likely that all factors have been 0'd out")
        message("Seek recourse.")
        return(which.min(params))
      }
        sd <- sd(bics)
        opt <- min(bics)
        in.range <-bics[(bics > (opt - sd)) & (bics < (opt+sd))]
        if(FALSE)
        {
          print("Optimal")
          print(opt)
          print("SD")
          print(sd)
          print(in.range)
          #top.indices <- which(bics %in% in.range)
          print("all BICs")
          print(bics)
          print("tops")
          print(top.indices)
          print(params)
          #optimal.l <- max(params[top.indices])
          print("Selecting params:")
          print(optimal.l)
        }
        top.indices <- which(bics %in% in.range)
        optimal.l <- max(params[top.indices])
        return(which(params == optimal.l))
    }
#Deals with cases if redundant scores.
#If these are at the upper extreme of the parameter list (likely occurs when all terms have been zeroed out), pick the SMALLEST parameter
#if these are at the lower extreme of the parameter list (likely occurs when the terms are fully dense), pick the largest parameter
  SelectBICFromIdenticalScores <- function(bic, params)
    {
      best.score <- min(bic)
      optimal.index <- which.min(bic)
      #message("Warning- BIC scores are unchanging for certain settings. This is likely evidence of no sparsity.")
     
      #Choose next lowest bic score index
      i <- sort(bic, index.return = TRUE)
      
      sorted.bic <- bic[i$ix]
      params.sorted.by.bic <- params[i$ix]
      matching.score.indices <- which(sorted.bic == best.score)
      
      #If the scores are on the  bigger end of scale
      if(all(min(params.sorted.by.bic[matching.score.indices]) > params.sorted.by.bic[-matching.score.indices]))
      {
        message("Suspect that scores are zeroing out the results, picking the smallest parameter with low BIC")
        optimal.param <- min(params.sorted.by.bic[matching.score.indices])
      } else if(all(max(params.sorted.by.bic[matching.score.indices]) < params.sorted.by.bic[-matching.score.indices]))
      {
        message("Suspect that scores are inducing no sparsity, picking the largest parameter with low BIC")
        optimal.param <- max(params.sorted.by.bic[matching.score.indices])
      } 
      else
      {
        #weird case, in the middle. #This means that all of them are equally goood?
        #I this case, I want to pick the most spare one actually
        print("Beware, unusual case...")
        #optimal.param <- params.sorted.by.bic[matching.score.indices][ceiling(length(matching.score.indices)/2)] #get the middle one
        optimal.param <- max(params.sorted.by.bic[matching.score.indices])
        print(sorted.bic)
        print(optimal.param)
        print(params.sorted.by.bic)
        print("")
        message("Unusual case with teh center, pick the middle. Likely swung too far above or below. Hope is lost :(")
        
      }
     
  which(params == optimal.param) #return the optimal index

    }
    
    #This gets the index for the optimal
  selectOptimalScoreIndex <- function(bic, params, oneSD, ndigits = 6)
  {
    bic <- round(bic, digits = ndigits)
    if(oneSD)#just testing our the one sd rule
    {
      optimal.index <- oneSDRule(bic,params)
    } else
    {
      best.score <- min(bic)
      optimal.index <- which.min(bic)
      if(length(which(bic == best.score)) > 1)
      {
        optimal.index <- SelectBICFromIdenticalScores(bic, params)
        if(length(optimal.index) > 1)
        {
          message("Warning: proposing redundant values.")
          optimal.index <- optimal.index[1]
        }
      }
    }
    optimal.index
  }
    #This function selects the optimal matrix and drops non-zero entries.
  selectOptimalInstance <- function(fit.data, bic, params, oneSD = FALSE)
  {
    
    optimal.index <- selectOptimalScoreIndex(bic, params, oneSD)
    optimal.matrix <- DropEmptyColumns(fit.data$fits[[optimal.index]][[1]])
    if(ncol(optimal.matrix) == 0)
    {
      optimal.matrix <- matrix(0, nrow = nrow(fit.data$fits[[optimal.index]][[1]]), ncol = 1 )
    }
    #SELECT NON-ZERO entries
  return(list("m" = optimal.matrix, "p" = params[[optimal.index]], "index" = optimal.index))
  }


    
getBICMatrices <- function(opath,option,X,W,all_ids, names, min.iter = 2, max.iter = 20, burn.in.iter = 4)
{
#If we get columns with NA at this stage, we want to reset and drop those columns at the beginning.
  burn.in.sparsity <- DefineSparsitySpaceInit(X, W, option, burn.in = burn.in.iter) #If this finds one with NA, cut them out here, and reset K; we want to check pve here too.

  option$K <- burn.in.sparsity$new.k
  consider.params <- SelectCoarseSparsityParams(burn.in.sparsity, burn.in.iter, n.points = 15)
  #things to record
  rec.dat <- list("alphas"=c(), "lambdas"=c(), "bic.a" = c(), "bic.l"=c(), "obj"=c(), 
                  "v.sparsity" = c(), "u.sparsity"=c(), "iter"=c(), "sd.sparsity.u" = c(), "sd.sparsity.v" = c(),
                  "alpha.s" = c(), "lambda.s" = c(), "Ks" = c(), "Vs" = list())
  #kick things off
  lambdas <- consider.params$lambdas
  alphas <- consider.params$alphas
	#option here to use the "burn in" sequence, but not surea bout that
  #optimal.v <- initV(X,W, option)
  optimal.v <- burn.in.sparsity$V_burn
  NOT.CONVERGED <- TRUE; i = 1
  #$Remove low PVE right now.
  optimal.v <- DropLowPVE(X,W,optimal.v, thresh = 0.01)
  while(i < max.iter & NOT.CONVERGED){
    print(i)
    #now fit U:
    rec.dat$Vs[[i]] <- optimal.v
	  #print(i)
    u.fits <- FitUs(X, W, optimal.v, alphas,option, weighted = TRUE)
    #record relevant info- always do 
    rec.dat$alphas <- c(rec.dat$alphas, alphas);  rec.dat$bic.a <- c(rec.dat$bic.a,u.fits$BIC) 
    #Pick the best choice from here, using threshold.
    optimal.iter.dat <- selectOptimalInstance(u.fits, u.fits$BIC, alphas)
    optimal.u <- optimal.iter.dat$m
    rec.dat$alpha.s <- c(rec.dat$alpha.s,optimal.iter.dat$p)

    rec.dat$U_sparsities = c(rec.dat$U_sparsities, matrixSparsity(optimal.u, ncol(X)));
    
    #now get the new lambdas for V:
    #Is this what we want?
    v.sparsity <- DefineSparsitySpace(X,W,as.matrix(optimal.u),"V", option)
    if(i == 1)
    {
      lambdas <- SelectCoarseSparsityParamsGlobal(v.sparsity, n.points = 15)
    }else
    {
      lambdas <- ProposeNewSparsityParams(v.fits$BIC, lambdas, v.sparsity, n.points = 7)
    }
    rec.dat$sd.sparsity.v <- c(rec.dat$sd.sparsity.v,sd(v.sparsity))
   
    #Now fit V!
    v.fits <- FitVs(X,W, optimal.u,lambdas,option, weighted = TRUE)
    #Pick the best choice from here, using threshold.
    optimal.iter.dat <- selectOptimalInstance(v.fits, v.fits$BIC, lambdas)
    optimal.v <- optimal.iter.dat$m
    #PercentVarEx(as.matrix(X)*as.matrix(W), v = optimal.v)
    #message("Updating K, iter ", i)
    rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v))
    rec.dat$lambda.s <- c(rec.dat$lambda.s,optimal.iter.dat$p)
    
    #Record new data
    rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
    rec.dat$lambdas <- c(rec.dat$lambdas, lambdas);  rec.dat$bic.l <- c(rec.dat$bic.l,v.fits$BIC) 

    #update the parameters for U based on the new V 
    u.sparsity <- DefineSparsitySpace(X,W,optimal.v, "U", option) #Here in case we hit the max.
    alphas <- ProposeNewSparsityParams(u.fits$BIC, alphas, (u.sparsity), n.points = 7)
    #Check convergence
    if(i > min.iter){
      
      NOT.CONVERGED <- !checkConvergenceBICSearch(i, rec.dat) #returns true if convergence is reached
    }
    if(ncol(optimal.v) == 1 & option$fixed_ubiq)
    {
      message("Only down to 1 factor, best be stopping now.")
      if(i < min.iter)
      {
        message("Zeroed-out the results very quickly, this suggests some kind of instability...")
      }
      NOT.CONVERGED <- FALSE
    }
    i = i+1
  }
  final.index <- i-1
  #TODO: add in drops for pve
  #This returns all the data from the last iteration
  #TODO: clean this up. This is very confusing.
  return(list("optimal.v" = optimal.v,"resid.var" = u.fits$resid_var,
              "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[final.index], "alpha"=rec.dat$alpha.s[final.index], "options" = option, 
              "K"= ncol(optimal.v), "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s))
}

#Convergence criteria for the BIC ssearch
#Converges when K is unchanging from one run to the next, and the percentage size change in the alpha/lambda paramters is less than 5%
#Might consider making this more generous- if it stays on the same log scale, then that is probably good enough....
checkConvergenceBICSearch <- function(index, record.data, conv.perc.thresh = 0.1)
{
  if(index > 5)
  {
    message("Late stage convergence, terminate soon....")
    #print(record.data$alpha.s)
    
    #print(record.data$lambda.s)
  }
  queries <- c(record.data$Ks[[index]] == record.data$Ks[[index-1]],
  abs(record.data$alpha.s[[index]] - record.data$alpha.s[[index-1]])/record.data$alpha.s[[index-1]] < conv.perc.thresh,
  abs(record.data$lambda.s[[index]] - record.data$lambda.s[[index-1]])/record.data$lambda.s[[index-1]] < conv.perc.thresh
  )
  return(all(queries))
}
#TODO: try with random, and with non-random.
gwasML_ALS_Routine <- function(opath, option, X, W, optimal.v, maxK=0)
{
#print(W)
  print(paste0("Starting at k:", option$K))
reg.run <- Update_FL(X, W, option, preV = optimal.v)
if(maxK != 0)
{
  if(ncol(reg.run$V)> maxK)
  {
    message("Parsing down to desired number of factors")
    message("More columns exist than specified....")
  }
  
  or <- sort(reg.run$PVE, decreasing =TRUE, index.return= TRUE)
  if(ncol(reg.run$V) < maxK)
  {
    message("Resulted in fewer than desired columns. Sorry.")
    maxK <- ncol(reg.run$V)
  }
  keep.cols <- or$ix[1:maxK]
  reg.run$V <- as.matrix(reg.run$V[,keep.cols])
  reg.run$U <- as.matrix(reg.run$U[,keep.cols])
  #print("Top PVEs:")
  #print(reg.run$PVE[keep.cols])
  #print("All PVEs")
  #print(reg.run$PVE)
  reg.run$PVE <- reg.run$PVE[keep.cols]
}

  save(reg.run, file = paste0(option$out,opath, "_gwasMF_iter.Rdata" ))
print(reg.run$V)

if(option$plots)
{
  o <- data.frame("rownames" = colnames(X), reg.run$V)
  write.table(o, file =  paste0(option$out,opath, ".factors.txt"), quote= FALSE, row.names = FALSE)
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
  title <- paste0("A=",reg.run$autofit_alpha[1], "Lambda=",reg.run$autofit_lambda[1]) 
  plotFactors(apply(reg.run$V, 2, function(x) x/norm(x, "2")), trait_names = o$rownames, title = title)
  ggsave(paste0(option$out,opath, ".factors.png"))
}

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
    message('Iniitializing X to the max -1')
    option$K <- ncol(X)-1
  }
  #Run the bic thing...
  bic.dat <- getBICMatrices(opath,option,X,W,all_ids, names)
  bic2.dat <- getBICMatrices(opath,option,X,W,all_ids, names)
  save(bic.dat, file = paste0(option$out,opath, "BIC_iter.Rdata" ))
  option <- bic.dat$options
  option$K <- bic.dat$K
  option$alpha1 <- bic.dat$alpha #Best from the previous. Seems to have converd..
  option$lambda1 <- bic.dat$lambda
  #I think the best thing to do is to initialize randomly, run it a few times, and either take the average or assess stability.
  #ret <- gwasML_ALS_Routine(opath, option, X, W, bic.dat$optimal.v)
  ret <- gwasML_ALS_Routine(opath, option, X, W, NULL) #randomly initialize here...
  ret
}

runStdPipeClean <- function(opath,args,alpha, lambda)
{
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
  option$alpha1 <- alpha
  option$lambda1 <- lambda
  #option$iter <- 5
  ret <- gwasML_ALS_Routine(opath, option, X, W, bic.dat$optimal.v)
  ret
}

gwasMFBIC <- function(X,W, snp.ids, trait.names, K=0, gwasmfiter =5)
{
  opath = "irrelevant"
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
  args <- defaultSettings(K=K)
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
  all_ids <-snp.ids; names <- trait.names
  initk <- option$K
  #1/16: testing out a new idea on the sims. Might be too much, but whatever.
  
  #if(initk == 0)
  #{
    message('Initializing X to the max -1')
    option$K <- ncol(X)-1
  #}
  #Run the bic thing...
  option$V <- FALSE
  #print("hard-coding in: K =5")
  #option$K = 5
  bic.dat <- getBICMatrices(opath,option,X,W,all_ids, names)
  #bic.dat2 <- getBICMatrices(opath,option,X,W,all_ids, names)
  #bic.dat3 <- getBICMatrices(opath,option,X,W,all_ids, names)
  #reportSimilarities(bic.dat, bic.dat2, bic.dat3)
  option <- bic.dat$options
  
 # option$K <- bic.dat$K
 option$K <-  ncol(X)-1
 #message("resetting to full k..")
  option$V <- TRUE
  option$alpha1 <- bic.dat$alpha #Best from the previous. Seems to have converd..
  option$lambda1 <- bic.dat$lambda
  #print("V preview:")
  #print(bic.dat$optimal.v)
  #NEED TO SEE IF GO RANDOM OR WHAT...
  #ret <- gwasML_ALS_Routine(opath, option, X, W, bic.dat$optimal.v, maxK = initk)
     #option$K <-  ncol(X)-1
  ret <- gwasML_ALS_Routine(opath, option, X, W, NULL, maxK=initk) #randomly initialize here...
  #ret2 <- gwasML_ALS_Routine(opath, option, X, W, NULL, maxK=initk) #randomly initialize here...
  #ret3 <- gwasML_ALS_Routine(opath, option, X, W, NULL, maxK=initk) #randomly initialize here...
  ret
}

defaultSettings <- function(K=0)
{
  args <- UdlerArgs()
  args$uncertainty <- ""
  args$gwas_effects <- ""
  args$nfactors <- K
  args$scale_n <- ""
  args$output <- "/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/matrix_simulations/RUN"
  opath <- "gwasMF"
  args$simulation <- TRUE
  args$converged_obj_change <- 0.05
  args$fixed_first <- TRUE #trying to see if this help- IT DOENS'T really appear to matter very much.
  args
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
  args$overview_plots <- FALSE
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


