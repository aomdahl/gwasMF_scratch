#Based on code from Yuan He, simplified



evaluate_error <- function(trueL, trueF, lp, fp){
  if(all(fp == 0) & all(lp == 0))
  {
    message("A zero case...")
    return(list("U_r2"=0, "V_r2"=0, "reordered_U"=lp,"reordered_V"= fp))
  }
  if(ncol(fp) < ncol(trueF)){
    dif = ncol(trueF) - ncol(fp)
    fp = cbind(fp, matrix(rep(0, nrow(fp) * dif), ncol = dif))
    lp = cbind(lp, matrix(rep(0, nrow(lp) * dif), ncol = dif))
  }
  rankK = ncol(trueF)
  suppressWarnings(library('combinat'))
  ordering = permn(seq(1,rankK))
  f_cor = rep(0, length(ordering))
  for(ord in seq(1, length(ordering))){
    f_cor[ord] = cor(as.vector(trueF), as.vector(fp[,ordering[[ord]]]))^2
  }
  
  l_cor = rep(0, length(ordering))
  for(ord in seq(1, length(ordering))){
    l_cor[ord] = cor(as.vector(trueL), as.vector(lp[,ordering[[ord]]]))^2
  }
  
  ord_sum = f_cor + l_cor
  ord = which.max(ord_sum)
  lp = lp[,ordering[[ord]]]
  fp = fp[,ordering[[ord]]]

  #lp = lp / matrix(rep(apply(lp, 2, function(x) max(abs(x))), nrow(lp)), nrow = nrow(lp), byrow = T)
  #fp = fp / matrix(rep(apply(fp, 2, function(x) max(abs(x))), nrow(fp)), nrow = nrow(fp), byrow = T)
  #colnames(fp) = seq(1,ncol(fp))
  #rownames(fp) = seq(1, nrow(fp))
  
  #fp[is.na(fp)] = 0 
  #lp[is.na(lp)] = 0
  
  l_r2 = cor(as.vector(trueL), as.vector(lp))^2
  f_r2 = cor(as.vector(trueF), as.vector(fp))^2
  

  return(list("U_r2"=l_r2, "V_r2"=f_r2, "reordered_U"=lp,"reordered_V"= fp))
}
