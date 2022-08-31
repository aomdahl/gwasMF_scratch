#New functions for new simulations:

##Whitening tools
estimate_null_correlation_simple <- function(z, z_thresh = 2, est_cor = TRUE) 
{
  #assumes z input
  max_absz = apply(abs(z), 1, max)
  nullish = which(max_absz < z_thresh)
  if (length(nullish) < ncol(z)) {
    stop("not enough null data to estimate null correlation")
  }
  nullish_z = z[nullish, ]
  Vhat = cov(nullish_z)
  if (est_cor) {
    Vhat = cov2cor(Vhat)
  }
  return(Vhat)
}

#From Guanghao fast asset
create_blocks <- function(cormat, cor_thr=0.2){
  # Hierarchical clustering
  corrdist <- as.dist(1-abs(cormat))
  hc <- hclust(corrdist)
  htree <- cutree(hc, h=1-cor_thr) # Criterion: corr>0.2
  block <- as.integer(names(table(htree))[table(htree)>=2])
  block <- lapply(block, function(x) names(htree)[htree==x])
  
  return(block)
}

blockify <- function(cormat, blocks)
{
  ret.mat <- diag(diag(cormat))
  colnames(ret.mat) <- colnames(cormat); rownames(ret.mat) <- rownames(cormat);
  if(length(blocks) != 0)
  {
    
    for(b in blocks)
    {
      if(length(b) >= 2)
      {
        for(i in b)
        {
          for(j in b)
          {
            ret.mat[i,j] <- cormat[i,j]
          }
        }
      }
    }
  }
  return(ret.mat)
  
}

#the whole thing.
blockCorrelationMatrix <- function(z,cor_thresh = 0.2, blocks = NULL)
{
  est <- estimate_null_correlation_simple(z) 
  if(any(is.na(blocks)))
  {
    blocks <- create_blocks(est, cor_thr=cor_thresh)
  }
 
  ret <- blockify(est, blocks)
  if(!isSymmetric(ret))
  {
    print("not symmetric?")
    print(ret)
  }
  return(ret)
}

#Functions to help process
#Get all GWAS
readInGWAS <- function(sim.results, dir)
{
  gwas.list <- list()
  for(i in 1:length(sim.results))
  {
    x <- sim.results[i]
    gwas.list[[i]] <- fread(paste0(dir,x))
  }
  return(gwas.list)
}

#Get PCA for all
pcaOnAll <- function(gwas.list, cols = "ALL")
{
  if(cols == "ALL")
  {
    lapply(gwas.list, function(x) svd(scale(x)))
  }else
  {
    lapply(gwas.list, function(x) svd(scale(x[,cols]))) #limit to just certain ones that contiain the data of interest.
  }

}

#Assess spearman correlation for all
  #spearmanVsPCs(true.null.pcs, c(0,0,0,1,1,1),trait.names = c("A1", "B1","C1", "A2", "B2", "C2") )
spearmanVsPCs <- function(run.name, pc.list, spearman.key,trait.names = c("A1", "B1", "A2", "B2"), stat.m = "spearman" )
{
  rho=list()
  p=list()
  for(i in 1:length(pc.list))
  {
    svs <- pc.list[[i]]
    pve <- svs$d^2/sum(svs$d^2)
    rownames(svs$v) <- trait.names
    #want to update this- looking at average
    sig.pves <- which(pve > mean(pve)) #Kaiser based method.
    #sig.pves <- which(pve > (1/length(svs$d)))
    rho[[i]] <- lapply(sig.pves, function(x) cor.test(y=spearman.key, x = svs$v[,x], method = stat.m)$estimate)
    p[[i]] <- lapply(sig.pves, function(x) cor.test(y=spearman.key, x = svs$v[,x], method = stat.m)$p.value)
  }
  ret <- NULL
  for(j in 1:length(rho))
  {
    for(i in 1:length(rho[[j]]))
    {
      ret <- rbind(ret, c(gsub(run.name[[j]],pattern = ".gwas.csv", replacement = ""), i, as.numeric(rho[[j]][[i]]), as.numeric(p[[j]][[i]])))
    }
  }
    colnames(ret) <- c("run", "PC", "rho", "pval")
    ret.frame <- data.frame(ret) %>% separate(run, sep = "\\.", into = c("seed", "rep", "overlap"))
    ret$overlap <- as.numeric(ret$overlap);ret$rho <- as.numeric(ret$rho);ret$pval <- as.numeric(ret$pval)
  return()
}

# Reconstruction error




# Pve change plot
pveChangePlot <- function(pc.list, list.names)
{
  null.pve <- do.call("rbind", lapply(1:length(pc.list), function(i) data.frame("pve" = pc.list[[i]]$d^2/sum(pc.list[[i]]$d^2)) %>% mutate("sim_id" = gsub(pattern=".gwas.csv", replacement = "",x = list.names[i])) %>% 
                                        tidyr::separate(sim_id, into=c("seed","sim_num", "overlap"), sep = "\\.", remove = TRUE) %>% mutate("pc"=1:length(pc.list[[i]]$d)))) %>% 
    mutate("perc_overlap" = round(as.numeric(overlap)/15000, digits = 2))
  ret.plot <- ggplot(null.pve, aes(x = as.factor(perc_overlap), y=pve )) + geom_boxplot() + 
    geom_smooth(method = "lm", se=FALSE, color="blue", aes(group="1")) + 
    facet_wrap(~pc,scales = "free_y" ) + theme_classic(15) + xlab("Percent overlap") + ylab("PVE")
  return(list("df" = null.pve, "plot" = ret.plot))
}

#Pairwise R2 plot:
#vp.null.sims
#look at 1-6
#lapply(1:6, function(x) cor(vp.null.zscores[[x]])[c])

pairwiseR2DF <- function(z.scores,sim.names, combs = matrix(c(c("A1","A1", "B1", "A2", "A2", "B2"),c("B1","C1", "C1", "B2", "C2", "C2")), ncol = 2))
{
  all.cors <- lapply(z.scores, function(x) cor(x)^2)
  combs.vect <- apply(combs, 1, function(x) paste0(x[1], ":", x[2]))
  tab.combs <- lapply(all.cors, function(c)  cbind(apply(combs, 1, function(x) paste0(x[1], ":", x[2])), c[combs]))
  joined.combs <- do.call("rbind", lapply(1:length(tab.combs), function(i) data.frame(tab.combs[[i]], "rep" = gsub(sim.names[i],pattern = ".gwas.csv", replacement = "")))) %>% 
    set_colnames(c("GWAS_entries", "R2", "run")) %>% 
    separate(run, sep = "\\.", into = c("seed", "rep", "overlap")) %>%
    mutate("cohort" = substr(GWAS_entries,start = 2, stop = 2) )
  joined.combs
}


#KS test tools
condKSTest <- function(dat, fcol, scol, thresh = 0.1)
{
  cond.tests <- dat[,scol][dat[,fcol] < thresh]
  ks.test(cond.tests, y="punif")
}

multiCondKS <- function(dat, pairs, pval = TRUE)
{
  if(pval) {return(lapply(pairs, function(x) condKSTest(dat, x[1], x[2])$p.value))}
  else
  {
    return(lapply(pairs, function(x) condKSTest(dat, x[1], x[2])$statistic))
  }
}


ksTestTab <- function(zin, sim.names, pairs = list(c("A1", "C1"), c("B1", "C1"), c("A2", "C2"), c("B2", "C2")), pair.names=c("A1:C1", "B1:C1", "A2:C2", "B2:C2"), pval = TRUE)
{
    ks.dat <- NULL
  for(tab in zin)
  {
    #make it p-values
    pvals <- data.frame(apply(abs(tab), 2, function(x) pnorm(-x)*2))
    if(pval)
    {
      std.unif <- apply(pvals, 2, function(x) (ks.test(x, y = "punif"))$p.value) 
    } else{
      std.unif <- apply(pvals, 2, function(x) (ks.test(x, y = "punif"))$statistic) 
    }
    
    conditionals <- unlist(multiCondKS(pvals, pairs, pval = pval))
    ks.dat <- rbind(ks.dat, c(std.unif, conditionals))
  }
  
  ks.df <- data.frame(ks.dat) %>% set_colnames(c(colnames(tab), pair.names)) %>% 
    mutate("run" = gsub(sim.names,pattern = ".gwas.csv",replacement = "")) %>% pivot_longer(cols = all_of(pair.names), "conditional_test") %>% separate(run, sep = "\\.", into = c("seed", "rep", "overlap"),remove = FALSE)
  ks.df$overlap <- as.numeric(ks.df$overlap)
  ks.df
}


#R2 performacne testing...
#pcas.each <- signal.simple.pcs
#sim.names <- signal.simple.sims
testFactorizationPerformance <- function(sim.names,pcas.each,npcs="detect", whitening = c())
{
  baseline <- which(grepl(sim.names, pattern = "\\.0.gwas.csv"))
  base.names <- sim.names[baseline]
  perf.tab <- NULL
  for(s in sim.names)
  {
    if(s %in% base.names)
    {
      message("Baseline run, skippping")
    }else
    {
      query <- str_split(s, pattern = "\\.")[[1]]
      base.match <- paste0(query[[1]],".", query[[2]], ".0.gwas.csv") #which is the correct one we want?
      base.pca <- pcas.each[[which(sim.names == base.match)]]
      query.pca <- pcas.each[[which(sim.names == s)]]
      if(npcs == "detect")
      {
        max.pc <- max(which(base.pca$d > mean(base.pca$d)))
      }else
      {
        max.pc = npcs
      }
      max.r2 <- evaluteFactorConstruction(as.matrix(base.pca$u[,1:max.pc]), as.matrix(base.pca$v[,1:max.pc]), as.matrix(query.pca$u[,1:max.pc]), as.matrix(query.pca$v[,1:max.pc]))
      if(length(whitening) > 0)
      {
        query.wht <- whitened.pcas.each[[which(sim.names == s)]]
        query.wht.raw <- whitened.pcas.raw[[which(sim.names == s)]]
        max.r2.wht <- evaluteFactorConstruction(as.matrix(base.pca$u[,1:max.pc]), as.matrix(base.pca$v[,1:max.pc]), as.matrix(query.wht$u[,1:max.pc]), as.matrix(query.wht$v[,1:max.pc]))
        max.r2.wht.raw <- evaluteFactorConstruction(as.matrix(base.pca$u[,1:max.pc]), as.matrix(base.pca$v[,1:max.pc]), as.matrix(query.wht.raw$u[,1:max.pc]), as.matrix(query.wht.raw$v[,1:max.pc]))
        perf.tab <- rbind(perf.tab, c(query[[1]], query[[2]], query[[3]], max.r2.wht, "pca_wht"))
        perf.tab <- rbind(perf.tab, c(query[[1]], query[[2]], query[[3]], max.r2.wht.raw, "pca_wht_raw"))
      }
      perf.tab <- rbind(perf.tab, c(query[[1]], query[[2]], query[[3]], max.r2, "pca_std"))

    }
  }
  data.frame(perf.tab) %>% set_colnames(c("seed", "rep", "overlap", "L_R2", "F_R2", "method")) %>% mutate("overlap" = as.numeric(overlap))
}
