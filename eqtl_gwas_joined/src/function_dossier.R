pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, Xmisc, cowplot, magrittr)
options(warn=-1)
###Visualization functions#####

#Cluster specification
clusterSpec <- function(joint.data.mat, ngwas, neqtl, dtCluster)
{
  if(dim(joint.data.mat)[1] > 1000000)
  {
    message('subsampling for clusters, data matrix too large...')
    clust.on <- joint.data.mat[sample(1:nrow(joint.data.mat), 10000),]
  } else
  {
    clust.on <- joint.data.mat
  }
  if(dtCluster)
  {
    #presumes eQTLs are first
    return(c(hclust(dist(t(clust.on[, 1:neqtl])))$order, hclust(dist(t(clust.on[, (neqtl + 1):(ngwas + neqtl)])))$order + neqtl))
    #c(hclust(dist(t(joint.data.mat[,1:ngwas])))$order, hclust(dist(t(joint.data.mat[, (ngwas + 1):(ngwas + neqtl)])))$order + ngwas)
  } else {
    return(hclust(dist(t(clust.on)))$order)
  }
}

#Correlation heatmap
#default cluster by data type to highlight groupings
#Do make it undirectional, simply pass in abs(joint.data.mat)
sepCorHeatmap <- function(joint.data.mat, all.names, ngwas = 55, neqtl = 49, dtCluster = TRUE)
{
  total_cor <- cor(joint.data.mat)
  #Make a nice heatmap
  #Order in a custom manner: by group, then hclust within groups...
  ordering <- clusterSpec(joint.data.mat, ngwas, neqtl, dtCluster)
  factors_nn <- data.frame(total_cor) %>% mutate("trait" = factor(all.names, levels = all.names[ordering]))
  colnames(factors_nn) <- c(all.names, "trait")
  nn <- tidyr::pivot_longer(factors_nn, cols = all_of(all.names), names_to = "trait2") %>%  arrange(value)
  factors_nn$src <- c(rep("green", neqtl), rep("black", ngwas))
  factors_nn <- factors_nn %>% arrange(trait)
  nn$trait2 <- factor(nn$trait2, levels =all.names[ordering])
  nn <- nn %>% arrange(trait, trait2)
  ggplot(nn, aes(trait, trait2, fill= value)) + geom_tile(color = "gray") +  
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab("Factors") + theme_minimal(10) +
    theme(axis.text.y = element_text(colour=factors_nn$src),axis.text.x = element_text(colour=factors_nn$src, angle = 45, hjust = 1)) + 
    scale_y_discrete(label=function(x) abbreviate(gsub(pattern = "_Analysis", replacement = "", x), minlength=25)) + 
    scale_x_discrete(label=function(x) abbreviate(gsub(pattern = "_Analysis", replacement = "", x), minlength=25))
}

#Make a PC plot
pcaPlotMaker <- function(L,joint.data.mat, all.names, title, ngwas = 55, neqtl = 49, dtCluster = TRUE){
  trait_names <- all.names
  new_names <- c(seq(1,ncol(L)), "trait")
  #Order in a custom manner: by group, then hclust within groups...
  ordering <- clusterSpec(joint.data.mat, ngwas, neqtl, dtCluster)
  #ordering <- orderFactors(F)
  factors_nn <- data.frame(L) %>% mutate("trait" = factor(trait_names, levels = trait_names[ordering]) )
  names(factors_nn) <- new_names
  nn <- tidyr::pivot_longer(factors_nn, cols = seq(1:ncol(L)), names_to = "x") %>%  arrange(value)
  nn$x <- as.factor(as.numeric(nn$x))
  factors_nn$src <- factors_nn$src <- c(rep("green", neqtl), rep("black", ngwas))
  factors_nn <- factors_nn %>% arrange(trait)
  ggplot(nn, aes(x, trait, fill= value)) + geom_tile(color = "gray") +  
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
    xlab("Factors") + theme_minimal(10) + theme(axis.text.y = element_text(colour=factors_nn$src)) + 
    ggtitle(title) + scale_y_discrete(label=function(x) abbreviate(x, minlength=30)) + 
    scale_x_discrete(guide = guide_axis(n.dodge=2))
  
  #theme(axis.text.x = element_text(vjust = grid::unit(c(-2, 2), "npc"))) 
}

#pcaBarPlotMaker(pca.sparse.gene,combined.mat.retry.numeric,full.name.list[-1], "Singular Vectors, joint PCA (sparse)", ngwas = 55, neqtl = 49)
pcaBarPlotMaker <- function(L,joint.data.mat, all.names, title, ngwas = 55, neqtl = 49, dtCluster = TRUE)
{
  trait_names <- all.names
  new_names <- c(seq(1,ncol(L)), "trait", "src")
  #Order in a custom manner: by group, then hclust within groups...
  ordering <- clusterSpec(joint.data.mat, ngwas, neqtl, dtCluster) #For large datasets, this step is slow. 
  coloration <- c(rep("eqtls", neqtl), rep("gwas", ngwas))
  #ordering <- orderFactors(F)
  
  factors_nn <- data.frame(L) %>% mutate("trait" = factor(trait_names, levels = trait_names[ordering]), "src" = coloration[ordering] )
  names(factors_nn) <- new_names
  nn <- tidyr::pivot_longer(factors_nn, cols = seq(1:ncol(L)), names_to = "x") %>%  arrange(value)
  nn$x <- as.factor(as.numeric(nn$x))
  nn$factors <- paste0("F", nn$x)
  
  p <- ggplot(nn, aes(x = trait, y = value, fill = src)) + geom_bar(stat='identity') + facet_wrap(~x) + 
    theme_minimal(10) + theme(axis.text.x=element_blank()) + xlab("Studies/traits") + ylab("Factor value")
  
  return(p)
}

#Analysis (PCA iterations)
customScale <- function(joint.data.mat, byCol = TRUE)
{
  if(!byCol) #scales across studies, by SNP
  {
    return(t(scale(t(joint.data.mat))))
  }
  return(scale(joint.data.mat) ) #I scale within a study, or by column
}

fullSpectralDecomp <- function(joint.data.mat, K = 50, byCol = TRUE)
{
  scaled <-  customScale(joint.data.mat, byCol)
  pca <- svd(x =scaled , nu = K, nv = K)
  
  png("test.png")
  dev.control(displaylist="enable")  
  plot(pca$d^2/sum(pca$d^2), col = "skyblue", pch = 19, xlab = "PC", ylab = "PVE", main = "PVE on joint PCA")
  abline(h=1/ncol(joint.data.mat), col="blue")
  p <- recordPlot()
  dev.off()
  
  list("svd" = pca, "scree" = p)
}

sparseSpectralDecomp <- function(joint.data.mat, K = 50, byCol = TRUE)
{
  library(ssvd)
  scaled <-  customScale(joint.data.mat, byCol)
  spca <- ssvd::ssvd(x = scaled, r = K)
  
  #for plotting
  png("test.png")
  dev.control(displaylist="enable") 
  plot(spca$d^2/sum(spca$d^2), col = "skyblue", pch = 19, xlab = "PC", ylab = "PVE", main = "PVE on Sparse Joint PCA")
  abline(h=1/ncol(joint.data.mat), col="blue")
  p <- recordPlot()
  dev.off()
  
  list("svd" = spca, "scree" = p)
}

#See which factors have overlappings study information across all PCs.
#returns a list indicating which SVs have overlap, all such studies overlapping, and a mapping between the two
mixtureCheck <- function(mat.in, all.names, top_n = 5, ngwas = 55, neqtl = 49)
{
  pcs <- c()
  overlap.names <- c()
  overlap <- list()
  counter = 1
  #The name picking is getting jacked up here for some reason...
  for(i in 1:ncol(mat.in))
  {
    ctr.score <- (mat.in[,i]^2)/sum(mat.in[,i]^2)
    over_avg <- which(ctr.score > (1/(ngwas + neqtl)))
    top <- intersect(order(-ctr.score)[1:top_n], over_avg)
    
    if(any(top <= neqtl) && any(top >= (neqtl + 1)))
    {
      print(paste0("PC number ", i, " has both GWAS and eqlt in the top ", top_n, " and above average."))
      pcs <- c(pcs, i)
      overlap.names <- unique(c(overlap.names, all.names[top]))
      overlap[[counter]] <- c(i, all.names[top])
      counter  = counter + 1
    }
  }
  return(list("pcs" = pcs, "studies" = overlap.names, "full" = overlap))
}

#Visualzie the mixture check object
#r is the mixturecheck object
mixtureViz <- function(r, all.names, spectral.decomp, num.col = 3, neqtl = 49, ngwas = 55 )
{
  #final viz
  message("Visualizing overlaps....")
  overlap.names <- r$studies
  #zoomin <- data.frame(spectral.decomp[(all.names %in% r$studies), r$pcs])
  #zoomin <- data.frame(spectral.decomp[, r$pcs])
  zoom.names <- all.names[(all.names %in% overlap.names)]
  trait.names <- all.names[(neqtl + 1):(neqtl+ngwas)]
  plotly <- list()
  if(length(r$pcs) > 0)
  {
    message("Should get overlaps for some factors....")
    for(j in 1:length(r$pcs))
    {
      i = r$pcs[j]
      zoomin <- spectral.decomp[, i]
      so <- order(zoomin^2,decreasing = TRUE)
      to.plot <- data.frame("vals" = ((zoomin[so]^2)/sum(zoomin^2))[1:10], "names" = all.names[so][1:10]) %>% 
        mutate("data_type" = ifelse(names %in% trait.names, "blue", "black"), "F" = "F1") %>% 
        arrange(vals^2)
      
      p <- ggplot(to.plot, aes(x = F, y = reorder(names,vals), fill = vals)) + 
        geom_tile(color = "gray" )+ theme_minimal(10) + 
        theme(axis.text.y = element_text(colour=to.plot$data_type), axis.text.x=element_blank()) + 
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(0,1))  + 
        xlab(paste0("SV", i))  + guides(fill=guide_legend(title=expression(v[ij]^2))) + 
        scale_y_discrete(label=function(x) abbreviate(x, minlength=30)) + ylab("") 
      legend <- get_legend(p)
      p <- p +  theme(legend.position = "none")
      #+ ylab("Study names")
      plotly[[j]] <- p
    }
    return(plot_grid(plotlist = plotly,nrow = ceiling(length(r$pcs)/num.col), ncol = num.col, legend))
  } else{return(ggplot() + theme_void())}
  
  
  
}


#run the full mixtures
#returns a mixtureCheck object and a plot of all the overlapping SVs
mixtureAnalysis <- function(spectral.decomp,all.names, ngwas, neqtl)
{
  r <- mixtureCheck(spectral.decomp, all.names, top_n = 5, ngwas = ngwas, neqtl =  neqtl)
  plot <- mixtureViz(r, all.names, spectral.decomp, neqtl = neqtl, ngwas = ngwas)
  message("Did you viz okay?")
  return(list("mixtures" = r, "plot" = plot))
}
##Correcting for covariates:
residualizeData <- function(adjust, joined.meta, joined.matrix, study.counts, all.names, args)
{
  if(adjust){
    message("going to do analysis via residuals, the straight up z scores.")
    #Here we regress by SNP, which will warimrW the effect of sample size per SNP, or one across all sNPs.
    #In order to the run the regression on the full matrix, we will subsample 10x and take the average.
    coef_mat <- matrix(NA, nrow = 10, ncol = 2)
    pval_mat <- matrix(NA, nrow = 10, ncol = 2)
    niter = 10
    nsamp = 1000
    for(i in 1:niter)
    {  
      counts.joint.data <- cbind(joined.meta$variant_id, joined.meta$gene_id, data.frame(joined.matrix)) %>% sample_n(nsamp) %>%
        pivot_longer(cols = colnames(.)[3:length(colnames(.))], names_to = "study" ) %>% 
        mutate("study" = cleanUpNames(study)) %>% left_join(., study.counts %>% mutate("study" = cleanUpNames(study)), by = "study")
      
      covars.correction <- lm(value ~ neff + src + 0, data = counts.joint.data) #no intercept term, effect sizes should be going through 0.
      
      coef_mat[i,] <- summary( covars.correction )$coef[,1]
      pval_mat[i,] <- summary( covars.correction )$coef[,4]
      counts.joint.data$resid <- covars.correction$residuals
      counts.joint.data$est <- covars.correction$fitted.values
      
    } #dropped intercept now
    
    keep_col <- which(colSums(pval_mat < 1) == 10)
    
    if(keep_col){
      betas <- as.matrix(colMeans(coef_mat)[keep_col])
      fitted <- as.matrix(all.study.counts %>% select(neff, src)) %*% (betas)
      
      #Helpful histogram of differences between original data and adjusted data (the double-sided plot was overkill)
      adj <-  matrix(fitted, nrow=nrow(joined.matrix), ncol=length(fitted), byrow=TRUE)
      joined.matrix = joined.matrix - adj
      png(paste0(args$output, "correction_hist.png"))
      hist(sample(adj, 1000), xlab = "Correction size", main = "Distribution of covariate corrections", 
           breaks = 30, col = "skyblue", border = "white")
      dev.off()
    }else{
      message("No significant effects associated with covariates; analysis will run as normal.")
    }
  } else{
    message("Not adjusting for covariates here.")
  }
  return(joined.matrix)
}

#########get study counts

getStudyCounts <- function(adjust, trait.names, eqtl.names)
{
  if(adjust)
  {
    message("Controlling for covariates of sample size (N-eff) and study type (gwas or eqtl)")
    gtex.count <- fread("/work-zfs/abattle4/ashton/reference_data/gtex_v8_sample_counts.txt", header = FALSE) #%>% mutate("tissue" = gsub(pattern = "-", replacement = "_", x = V1))
    gwas.count <- fread("/work-zfs/abattle4/ashton/snp_networks/scratch/eqtl_gwas_joined/fourth_pass_universal/gwas_extracts/herit_list.txt")
    trait.counts <- left_join(trait.names %>% rename("phenotype" = V1), gwas.count, by = "phenotype") %>% select(V2, Neff) %>% mutate("src" = 1)  %>% magrittr::set_colnames(c("study", "neff", "src")) #lots more information here, but for now keep it simple
    tiss.counts <- left_join(data.frame("V1" = eqtl.names), gtex.count, by = "V1") %>% mutate("src" = 0) %>% magrittr::set_colnames(c("study", "neff", "src"))
    #0 is eqlt, 1 is gwas
    all.study.counts <- rbind(trait.counts, tiss.counts)
  } else{all.study.counts = NULL}
  return(all.study.counts)
}
########General i/o
speedRead <- function(path)
{
  path_ = gsub(pattern = "rand_", replacement = "rand.", x = path)
  fread(path) %>% mutate("src" = str_split(string= basename(path_), pattern = "\\.")[[1]][1])
}



maxfValue <- function(dat) {
  #Note we need to go 2:50 since it drops the grouping variable
  data.frame(t(apply(dat[,2:ncol(dat)], 2, function(x) x[which.max(abs(x))]))) #get the max snp-gene pair, snp fixed.
}



s#Options are "RANDOM", "MAX_AVG_Z", "CENTER_Z", "MAX_GENE"
#parseEQTLs(args$eqtl_source, args$eqtl_selection)
parseEQTLs <- function(dir, sel.criteria, num = "ALL")
{
  #Parse in the eQTL extracted data
  #dir <- "/work-zfs/abattle4/ashton/snp_networks/scratch/eqtl_gwas_joined/fourth_pass_universal/eqtl_extracts/500kb_0.1r2/"
  files <- paste0(dir,"/", list.files(dir, pattern = "*.extract.tsv"))
  if(num != "ALL")
  {
    files <- files[1:num]
  }
  
  #Default just get z_score
  eqtl.long <- do.call("bind_rows", lapply(files, speedRead)) %>% 
    magrittr::set_colnames(c("gene_id", "variant_id", 
                             "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se", "src")) 
  
  #If really really long....
  #Subset all the SNPs to the minimum length one
  #counts <- eqtl.long %>% group_by(src) %>% dplyr::summarize("c" = n()) %>% slice(which.min(c))
  #src_ = unlist(counts[1,1])
  #filt.list <- (eqtl.long %>% filter(src == src_))$variant_id
  #t <- eqtl.long %>% filter(variant_id %in% filt.list)
  message("Pivoting wide")
  #use tidyfast to do this faster..
  #faster with dt_pivot_wider()
  eqtl.zscores <- data.table(eqtl.long %>%  mutate("z_score" = slope/slope_se) %>% 
    select(gene_id, variant_id, z_score, src) %>% 
    pivot_wider(id_cols = c("gene_id", "variant_id"), names_from = "src", values_from = z_score))
  # This should technically not be necesary. At any rate, our current approach is too slow%>%     stdGenomicFilter(tab = ., col_name = "variant_id")
  
  #Now, specify the filtering mechanism
  if(sel.criteria == "MAX_GENE")
  {
    message("For convenience, we replace missing entries with 0 in the eQTL data.")
    eqtl.zscores[is.na(eqtl.zscores)] <- 0
    message("Also adding a blank column of genes, since its variable")
  }
  selected.eqtl <- switch(  
    sel.criteria,  
    "RANDOM"= eqtl.zscores %>% group_by(variant_id) %>%  slice_sample(n = 1),
    "MAX_AVG_Z"= eqtl.zscores %>% mutate("avg" = rowMeans(abs(.[3:ncol(.)]),na.rm = TRUE)) %>% group_by(variant_id) %>%  
      slice(which.max(avg)) %>% select(-avg),
    "MEDIAN_Z"= eqtl.zscores %>% mutate("avg" = rowMeans(abs(.[3:ncol(.)]), na.rm = TRUE)) %>% group_by(variant_id) %>%  
      slice(which.min(abs(avg-median(avg)))) %>% select(-avg),
    "MAX_GENE" = eqtl.zscores %>% group_by(variant_id) %>% group_modify(~ maxfValue(.x)) #%>% mutate("gene_id" = "var")
  )  %>% ungroup() %>% arrange(variant_id)
  
  return(selected.eqtl)
}

if(FALSE)
{
  
  #code to do it
  topByGene <- function(dat, var, gene){
    dat[is.na(dat)] = 0 #this is important, ensures we get some gene picked for each one.
    sizes <- list()
    ngenes.to.choose.from <- c()
    iter = 1
    extract.query <- NULL
    all.snps <- unique(unlist(dat[var]))
    #alt version
    toy <- eqtl.zscores[1:1000,]
    #just get the max value out
    
    test <- toy %>% group_by(variant_id) %>% group_modify(~ maxfValue(.x))
    #Get the max gene out
    maxfGene <- function(dat) {
      gene = "gene_id"
      #print(colnames(dat))
      #Note we need to go 2:50 since it drops the grouping variable
      sel_out <- unlist(apply(dat[,2:50], 2, function(x) which.max(abs(x)))) #get the max snp-gene pair, snp fixed.
      unlist(dat[gene])[(sel_out)]
    }
    #Important note: if our "group_by" value is a string, the result is sorted
    #Check now to ensure we actually get the max...
    aa <- (toy %>% filter(variant_id == "1:998726:A:G"))$Brain_Amygdala
    aa[which.max(abs(aa))]
    #Checks out.
    eqtl.zscores[is.na(eqtl.zscores)] <- 0
    full <- eqtl.zscores %>% group_by(variant_id) %>% group_modify(~ maxfValue(.x))
    save(full, file = "/work-zfs/abattle4/ashton/snp_networks/scratch/eqtl_gwas_joined/fourth_pass_universal/snp_fixed/500kb_0.1r2.eqtl.RData")
    #Next step: build a lookup key
    sorted.snps <- sort(unique(toy$variant_id))
    t <- do.call("rbind", lapply(1:length((test)), function(x) data.frame("genes" = as.character(test[[x]]), "var" = as.character(sorted.snps[x]))))
    
    for(snp in all.snps)
    {
      #faster to group by snp, arrange by max gene, select top?
      curr <- dat[dat[var] == snp,]
      ngenes.to.choose.from <- c(ngenes.to.choose.from, nrow(curr))
      sel_out <- apply(curr[,3:51], 2, function(x) which.max(abs(x))) #get the max snp-gene pair, snp fixed.
      query.genes <- unlist(curr[gene])[(sel_out)]
      #set up snp-gene - tissue pairs to extract for each thing...
      extract.query <- rbind(extract.query, data.frame("SNP_GENE" = paste0(snp,"_", query.genes), "Tissue" = 2:50 ))
      iter = iter + 1
      progress(iter, max.value = length(all.snps))
    }
    head(extract.query)
  }
}

handleNAs <- function(tab, method)
{
  message(paste0("NAs will be handled with ", method, " method."))
  ret <- switch(  
    method,  
    "ZERO"= {tab[which(is.na(tab))] <- as.double(0); tab},
    "DROP"= tab %>% drop_na(), #not implemented- need to track which indices get dropped.
    "SNP_MEAN_IMPUTE"= tab #not yet implemented
  )  
  return(ret)
}

joinData <- function(eqtls.z, gwas.z, method = "ZERO")
{
  #join the data
  joint.data <- left_join(eqtls.z, gwas.z, by = "variant_id") %>% arrange(hg38) 
  joint.meta.data <- joint.data %>% select(hg38, hg19, variant_id, ref, alt, gene_id)
  joint.data <- as.matrix(joint.data %>% select(-hg38, -hg19, -variant_id,-ref,-alt, -gene_id)) %>% handleNAs(., method)
  ret <- list("mat" = joint.data, "meta" = joint.meta.data)
  return(ret)
}

#helper to save all the images we generate.
alignGWASNames <- function(trait.tab, local)
{
  return((data.table("V1" = local) %>% left_join(., trait.tab) %>% filter(V1 %in% trait.tab$V1))$V2)
}

savePlots <- function(names, plot.objs, main.path)
{
  for(i in 1:length(names))
  {
    if(grepl("PLOT",names[i]))
    {
      png(file=paste0(main.path, names[i], ".png"))
      print(plot.objs[[i]])
      dev.off()
    }else{
      ggsave(plot = plot.objs[[i]], filename = paste0(main.path, names[i], ".png"),width = 9, height = 9)
    }
  }
}

#basic tool to filter according
#@param tab: the table
#@param col_name: the name of teh column containing id, here chr:pos:ref:alt
stdGenomicFilter <- function(tab, col_name)
{
  #Does the following in this order:
  # Drop duplicates (both copies)
  # Remove indels
  # Remove ambigs
  ambig <- c("AT", "TA", "GC", "CG")
  ret <- tab %>% separate(!!sym(col_name), remove = FALSE, into = c("chr", "pos", "ref", "alt")) %>% 
    mutate("reg_id" = paste0(chr, ":", pos)) %>% group_by(reg_id, gene_id) %>% filter(n() == 1) %>% ungroup() %>%
    mutate("joint" = paste0(ref,alt)) %>% filter(nchar(joint) == 2) %>%
    filter(!(joint %in% ambig)) %>% select(-chr, -pos, -ref, -alt, -joint, -reg_id)
  return(ret)
}

cleanUpNames <- function(xin)
{
  gsub(x = xin, pattern = "\\)", replacement = "") %>% 
    gsub(x = ., pattern = "\\(", replacement = "") %>%
    gsub(x = ., pattern = "\\.", replacement = "") %>%
    gsub(x = ., pattern = "\\:", replacement = "") %>%
    gsub(x = ., pattern = "X(\\d)", replacement = "\\1") %>%
    gsub(x = ., pattern = "\\,", replacement = "") %>%
    gsub(x = ., pattern = "-", replacement = "")
  #check
}
#fullFactorizationAnalysis(joint.data$mat,all.names, ngwas, neqtl, "signed")
#make sure to include "PLOT" in the r std plot format, makessure we wirte out correctly later.
fullFactorizationAnalysis <- function(data.mat,all.names, ngwas, neqtl, h_handle)
{
  pca <- fullSpectralDecomp(data.mat,byCol = TRUE)
  x <- scale(data.mat)
  #top.pcs <- pca$svd$u[,1:2] %*% (t(pca$svd$u[,1:2]) %*% x)
  #Also removing the first 2 PCs, and trying the heatmap on that.
  top.pcs.straight <- pca$svd$u[,1:2] %*% diag(pca$svd$d[1:2]) %*% t(pca$svd$v[,1:2])
  residualized = x - top.pcs.straight
  
  heatmap <- sepCorHeatmap(data.mat, all.names, ngwas = ngwas, neqtl = neqtl, dtCluster = TRUE)
  heatmap.res <- sepCorHeatmap(residualized, all.names, ngwas = ngwas, neqtl = neqtl, dtCluster = TRUE)
  
  #For fun- try removing  a few more PCs...
  top.pcs.straight <- pca$svd$u[,1:4] %*% diag(pca$svd$d[1:4]) %*% t(pca$svd$v[,1:4])
  residualized = x - top.pcs.straight
  heatmap.res.4 <- sepCorHeatmap(residualized, all.names, ngwas = ngwas, neqtl = neqtl, dtCluster = TRUE)
  
  #Vanilla PCA
  pca <- fullSpectralDecomp(data.mat,byCol = TRUE)
  pca.heatmap.dense <- pcaPlotMaker(pca$svd$v,data.mat, all.names, "Singular Vectors, joint PCA (std)", ngwas = ngwas, neqtl = neqtl)
  pca.barplot.dense <- pcaBarPlotMaker(pca$svd$v,data.mat, all.names, "Singular Vectors, joint PCA (std)", ngwas = ngwas, neqtl = neqtl)
  
  #Sparse PCA
  pca.sparse <- sparseSpectralDecomp(data.mat,byCol = TRUE)
  pca.heatmap.sparse <- pcaPlotMaker(pca.sparse$svd$v, data.mat, all.names, "Singular Vectors, joint PCA (sparse)", ngwas = ngwas, neqtl = neqtl)
  pca.barplot.sparse <- pcaBarPlotMaker(pca.sparse$svd$v,data.mat, all.names, "Singular Vectors, joint PCA (sparse)", ngwas = ngwas, neqtl = neqtl)
  
  #mixture analysis
  message("Beginning PCA mixture analysis...")
  mix <- mixtureAnalysis(pca$svd$v, all.names, ngwas = ngwas, neqtl = neqtl)
  message("Beginning sparse PCA mixture analysis...")
  mix.sparse <- mixtureAnalysis(pca.sparse$svd$v, all.names, ngwas = ngwas, neqtl = neqtl)
  ret <- list()
  ret$plots <- list(heatmap, pca.heatmap.dense, pca.barplot.dense, pca$scree, pca.heatmap.sparse, 
                    pca.barplot.sparse, pca.sparse$scree, mix$plot, mix.sparse$plot, heatmap.res,heatmap.res.4)
  ret$fnames <- paste0(c("corr_heatmap", "SVs_dense","barplot_dense", "scree_PLOT_dense", "SVs_sparse",
                         "barplot_sparse","scree_PLOT_sparse", "mixture_dense", "mixture_sparse", "corr_resid2_heatmap","corr_resid4_heatmap"), "_", h_handle)
  
  return(ret)
}

#Speciify which PCs to remove
#abs(joint.data.mat),all.names, c(1), 55, 49, "unsigned_custm"
customHeatmapAnalysis <- function(data.mat,all.names, droppers, ngwas, neqtl, h_handle)
{
  pca <- fullSpectralDecomp(data.mat, byCol = TRUE)
  x <- scale(data.mat)
  #top.pcs <- pca$svd$u[,1:2] %*% (t(pca$svd$u[,1:2]) %*% x)
  #Also removing the first 2 PCs, and trying the heatmap on that.
  if(length(droppers) == 1)
  {
    diag = pca$svd$d[droppers]
    top.pcs.straight = diag * pca$svd$u[,droppers] %*% t(pca$svd$v[,droppers])
  }else
  {
    top.pcs.straight <- pca$svd$u[,droppers] %*% diag(pca$svd$d[droppers]) %*% t(pca$svd$v[,droppers])
  }
  
  residualized = x - top.pcs.straight
  
  heatmap <- sepCorHeatmap(data.mat, all.names, ngwas = ngwas, neqtl = neqtl, dtCluster = TRUE)
  heatmap.res <- sepCorHeatmap(residualized, all.names, ngwas = ngwas, neqtl = neqtl, dtCluster = TRUE)
  
  ret <- list()
  ret$plots <- list(heatmap,heatmap.res)
  ret$fnames <- paste0(c("corr_heatmap", "corr_resid_custom_heatmap"), "_", h_handle)
  ret$res <- residualized
  
  return(ret)
}

