#Manual check of results.....
library("optparse")
library("data.table")
library("magrittr")
library("dplyr")
library("combinat")
source("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/src/evaluateSimR2.R")
library(NatParksPalettes)
option_list <- list( 
  make_option(c("-o", "--output"), default = '', type = "character", 
              help="Output dir + handle"),
  make_option(c("-s", "--sim_path"), default = '', type = "character", 
              help="File containing all of the paths of simulations you wish to go into this plot."),
  make_option(c("-p", "--plot"), default = FALSE, type = "logical", action = "store_true", 
              help="specify this if you want to plot."),
  make_option(c("-y", "--yaml"), default = '', type = "character", 
              help="Path to the settings file.")
)

#debug interactive
t = c("--output=scratch16-abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/udler_2_k4_max_all//factorization_results/summary",
      "--sim_path=scratch16-abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/udler_2_k4_max_all//factorization_results/",
      "--yaml=scratch16-abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/udler_2_k4_max_all/udler_2_k4_max_all.yml")
args <- parse_args(OptionParser(option_list=option_list), args = t)



yml <- read.table(args$yaml, sep = ",")
methods.run <- ((yml %>% filter(V1 == "test_methods"))$V2 %>% strsplit(.,":" ))[[1]]
true.loadings <- as.matrix(fread((yml %>% filter(V1 == "loadings"))$V2))
true.factors <-  as.matrix(fread((yml %>% filter(V1 == "factors"))$V2))
niter <- as.numeric((yml %>% filter(V1 == "iter"))$V2)
s <- args$sim_path
m <- "PCA"
#m<-"sPCA"
i = 1
pred.factors <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".factors.txt")))
pred.loadings <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".loadings.txt")))
reconstruction <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".recon_error.txt")))

#This actually seems like a great way to do it. Maybe should be doing this instead.

ordering = permn(4)
library(gtools)
sign.ordering <- gtools::permutations(2, 4, c(1,-1), repeats.allowed = TRUE)


getallR2 <- function(sign.ordering, ordering, pred, actual, reverse = FALSE)
{
  nreps = nrow(actual)
  all.options.factors <- apply(sign.ordering, 1, function(s) lapply(ordering, function(x) pred[,x] * do.call("rbind", lapply(1:nreps, function(x) s))))
  res.mat <- matrix(NA,nrow(sign.ordering),length(ordering))
  #For each of those, estimate the R2
  for(signopt in 1:nrow(sign.ordering))
  {
    res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) cor(as.vector(x), as.vector(actual))^2))
    if(reverse)
    {
      res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) cor(as.vector(actual),as.vector(x))^2))
    }
  }
  return(res.mat)
}



#on the F side....
#run both iterating true and pred.
#create every combination of sign and order.
f.scores <- getallR2(sign.ordering, ordering, pred.factors, true.factors)
#On the l side
l.scores <- getallR2(sign.ordering, ordering, pred.loadings, true.loadings)
tot = l.scores + f.scores

l.scores[which.max(tot)]
f.scores[which.max(tot)]

#Do these match the outputted?
#This is actually a way cleaner way to accomplish the same thing...

#try in reverse...

f.scores.b <- getallR2(sign.ordering, ordering, as.matrix(true.factors),pred.factors)
#On the l side
l.scores.b <- getallR2(sign.ordering, ordering, as.matrix(true.loadings),pred.loadings)
tot.b = l.scores.b + f.scores.b

max(tot.b)
max(tot)


#Try with the entries swapped...
f.scores.b.r <- getallR2(sign.ordering, ordering, as.matrix(true.factors),pred.factors,reverse = TRUE)
#On the l side
l.scores.b.r <- getallR2(sign.ordering, ordering, as.matrix(true.loadings),pred.loadings, reverse = TRUE)
tot.b.r = f.scores.b.r+l.scores.b.r
#The global max ought to be the same as the global max elsewhere- at last the same combos.
which.max(tot.b.r)
max(tot.b.r)
max(tot.b)
max(tot)
#22, sign 6, order 2
which.max(tot) #23, sign 7, order 2

#Bigger difference in F, look there.
f.scores.b[22]
f.scores[23]

sign.ordering[6,]
sign.ordering[7,]
ordering[[2]]
true.factors[,ordering[[2]]]
pred.factors

pred.factors[,ordering[[2]]]
true.factors
signExpand <- function(t, N)
{
  do.call("rbind", lapply(1:N, function(x) t))
}

cor(as.vector(pred.factors[,ordering[[2]]] * signExpand(sign.ordering[7,], 10)) , as.vector(true.factors))^2

cor(as.vector(true.factors[,ordering[[2]]] * signExpand(sign.ordering[6,], 10)) , as.vector(pred.factors))^2

message("original")
cor((pred.factors[,ordering[[2]]] * signExpand(sign.ordering[7,], 10)) , (true.factors))^2

message("backwards")
cor((true.factors[,ordering[[2]]] * signExpand(sign.ordering[6,], 10)) , (pred.factors))^2

#Everything looks right and seems to check out. Not really sure what is going on here.

#The errors are close enough that I don't think it matters, but might be worth checking some more

#Empirically, it appears that while differences do exist, they are so very small

#leaving a big question mark on this, moving on.

#check the correlation matrix C over time...