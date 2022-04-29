#Generate phenotypes

#snp.weights- which snps are active?
#build your phenotypes using specified snps and genotypes
phenotypeBuilder <- function(snp.weights, genotype, noise, type = "1+2=3")
{
  ret.pheno = list()
  if(type == "1+2=3")
  {
    #If 3 is a linear combination of 1 and 2
    #Model 1 and 2 as independent
    ret.pheno[[1]] = genotype %*% snp.weights[[1]] + rnorm(nrow(genotype), 0, sd = noise)
    ret.pheno[[2]] = genotype %*% snp.weights[[2]] + rnorm(nrow(genotype), 0, sd = noise)
    ret.pheno[[3]] = 0.5*ret.pheno[[1]] +  0.5*ret.pheno[[2]] + rnorm(nrow(genotype), 0 ,sd = noise)
  }
  if(type == "1,2,3") #allows for some overwlapping snps between them.
  {
    #If 3 is a linear combination of 1 and 2
    #Model 1 and 2 as independent
    ret.pheno[[1]] = genotype %*% snp.weights[[1]] + rnorm(nrow(genotype), 0 ,sd = noise)
    ret.pheno[[2]] = genotype %*% snp.weights[[2]] + rnorm(nrow(genotype), 0 ,sd = noise)
    ret.pheno[[3]] = genotype %*% snp.weights[[3]] + rnorm(nrow(genotype), 0 ,sd = noise)
  }
  if(type == "1,2") #allows for some overwlapping snps between them.
  {
    #If 3 is a linear combination of 1 and 2
    #Model 1 and 2 as independent
    ret.pheno[[1]] = genotype %*% snp.weights[[1]] + rnorm(nrow(genotype), 0 ,sd = noise)
    ret.pheno[[2]] = genotype %*% snp.weights[[2]] + rnorm(nrow(genotype), 0 ,sd = noise)
  }
  if(substr(type,1,1) == "0")
  {
    s = str_split(type,pattern = ",")[[1]]
    for(i in 1:length(s))
    {
      if(s[i] == "0")
      {
        ret.pheno[[i]] = rnorm(nrow(genotype),0,noise)
      }else
      {
        message("not yet implemented")
      }
      
    }
  }
  return(ret.pheno)
}



relatedOverlapPheno <- function(n.snps, n.groups, N.tot = 10000, genotypes=genotypes.tab.2, 
                                pheno = "1+2=3", size = 0.5, active_snps = 30, overlap_snps = 25)
{
  print(N.tot)
  genotypes <- genotypes[,1:n.snps]
  #message(interval_size)
  nonoverlapping <- active_snps-overlap_snps
  active.snps <- list()
  if(pheno =="1+2=3" )
  {
    
    active.snps[[1]] <- c(rep(size, 20), rep(0, n.snps-20))
    active.snps[[2]] <- c(rep(0, 20), rep(size, 20), rep(0, n.snps-40))
    active.snps[[3]] <- c(rep(0, n.snps))
  } else if(pheno == "1,2,3") #case with overlap- 2 similar traits
  {
    active.snps[[1]] <- c(rep(size, 40), rep(0, n.snps-40))
    active.snps[[2]] <- c(rep(0, 20), rep(size, 40), rep(0, n.snps-60))
    active.snps[[3]] <- c(rep(0, n.snps))
  } else if(pheno == "1,2") #case with overlap- 2 similar traits
  {
    active.snps[[1]] <- c(rep(size, active_snps), rep(0, n.snps-active_snps))
    active.snps[[2]] <- c(rep(0, nonoverlapping), rep(size, active_snps), rep(0, n.snps-(active_snps + nonoverlapping)))
  } else if(substr(pheno,1,1) == "0") #first element is null
  {
    s = str_split(type,pattern = ",")[[1]]
    for(i in 1:length(s))
    {
      if(s[i] == "0")
      {
        active.snps[[i]] <- c(rep(0, n.snps))
      }else
      {
        message("not yet implemented")
      }
      
    }
  }
  else{
    message("forgot to specify type...")
  }
  phenotypeBuilder(active.snps, genotypes, 0.5,  type = pheno) #this returns a list of all phenotypes for all individuals.
  
}

relatedOverlapBETTER <- function(n.snps, n.overlap, n.groups, y, N.tot = 10000, 
                           genotypes=genotypes.tab.2, fixed.count = 5000, pheno = "1+2=3")
{
  if(n.overlap > 0)
  {
    No <- 1:n.overlap
  }else
  {
    No <- NULL
  }
  interval_size = floor((N.tot-n.overlap)/n.groups)
  g <- list()
  t = n.overlap
  
  if(t >= n.snps)
  {
    t = n.overlap
  }
    #the number of phenotypes must correspond with the number of traits.
    if(fixed.count  == 0)
    {
      for(i in 1:n.groups)
      {
        message(paste0("group " , i))
        if(interval_size > 0)
        {
          cohort_length <- c(No, (t+1):(t+interval_size))
          print(length(cohort_length))
        } else
        {
          cohort_length  <- No
        }
        g[[i]] <- sapply(1:n.snps, function(j) summary(lm(y[[i]][cohort_length] ~ genotypes[cohort_length, j]))$coef[6])
        t= t+interval_size
      }
    }else #fixed count implies that each cohort has the exact same number of indivudals
    {
      No <- 1:n.overlap
      add <- fixed.count - n.overlap
      for(i in 1:n.groups)
      {
        start_point = (i-1) * fixed.count -((i>2) * n.overlap)
        cohort_length <- c(No, (start_point + ifelse(i == 1, n.overlap, 0)):(start_point + ifelse(i == 1, n.overlap, 0) + add))
        if(cohort_length[fixed.count] > N.tot)
        {
          print("Error")
        }
        print(length(cohort_length))
        g[[i]] <- sapply(1:n.snps, function(j) summary(lm(y[[i]][cohort_length] ~ genotypes[cohort_length, j]))$coef[6])
      }
    }
    
    return(do.call("cbind", g))
    }

