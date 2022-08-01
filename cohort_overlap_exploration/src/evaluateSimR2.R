#### Function to evaluate error between predicted L and F and true L and F
#### Taken from Yuan He's code with some modification.

#Changes:
  #Sig_hits not considered 

evaluteFactorConstruction <- function(trueL, trueF, lp, fp){
  #fill in missing columns with 0s
  if(ncol(fp) < ncol(trueF)){
    message("missing factors, will set to 0...")
    dif = ncol(trueF) - ncol(fp)
    fp = cbind(fp, matrix(rep(0, nrow(fp) * dif), ncol = dif))
    lp = cbind(lp, matrix(rep(0, nrow(lp) * dif), ncol = dif))
  }
  rankK = ncol(trueF)
  suppressWarnings(library('combinat'))
  if(K > 9)
  {
    message("There are too many latent factors for comprehensive examination of options. Do a greedy search.")
    return()
  }else
  {
    #get all possible permutations in F and L
    ordering = permn(rankK)
    f_cor = rep(0, length(ordering))
    for(ord in seq(1, length(ordering))){
      f_cor[ord] = cor(as.vector(trueF), as.vector(fp[,ordering[[ord]]]))
    }
    
    l_cor = rep(0, length(ordering))
    for(ord in seq(1, length(ordering))){
      l_cor[ord] = cor(as.vector(trueL), as.vector(lp[,ordering[[ord]]]))
    }
    
    #Look at the sum of correlations overall (taking absolute value)
    #Note that there may be an exceptional case, where they have opposite directions of correlation, and this is missed
    ord_sum = (f_cor + l_cor)^2 #the squared sum, becasue effects may be positive or negative
  
    ord = which.max(ord_sum)  #which combination is highest overall
    #order columsn by this
    lp = lp[,ordering[[ord]]]
    fp = fp[,ordering[[ord]]]
    
    #Omitting these steps- not sure about their purpose, nor their usefulness here where we aren't semi-non-negative.
    #lp = lp / matrix(rep(apply(lp, 2, function(x) max(abs(x))), nrow(lp)), nrow = nrow(lp), byrow = T) #For some reason, we are scaling by the max value per row in each entry. I am not really sure why.
    #fp = fp / matrix(rep(apply(fp, 2, function(x) max(abs(x))), nrow(fp)), nrow = nrow(fp), byrow = T)
    colnames(fp) = seq(1,ncol(fp))
    rownames(fp) = seq(1, nrow(fp))
    
    fp[is.na(fp)] = 0 
    lp[is.na(lp)] = 0
    
    l_corr = cor(as.vector(trueL), as.vector(lp))
    f_corr = cor(as.vector(trueF), as.vector(fp))

    return(c(l_corr, f_corr))
  }
  

}
