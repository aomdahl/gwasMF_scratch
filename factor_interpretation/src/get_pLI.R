#INTERROGATING PLI SCORES
#DONWLOADED FROM https://gnomad.broadinstitute.org/downloads#v4-constraint
#WITH wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv

#Its in :  /data/abattle4/aomdahl1/reference_data/gnomad.v4.1.constraint_metrics.tsv


#pli_scores <- fread("/data/abattle4/aomdahl1/reference_data/gnomad/gnomad.v4.1.constraint_metrics.tsv")

#can do this acrtoss multiple factors

getConstraintDist <- function(input.genes, pli_scores = NULL)
{
  if(is.null(pli_scores))
  {
    pli_scores <- getPLIDat()
  }
  #This eliminates redundant matches- we want to account for redundancy, so use a left join
  input.look <- data.frame("gene" = input.genes)
  #matches <- (pli_scores %>% filter(gene %in% input.genes))
  #Need to drop ones where no score is available
  matches <- left_join(input.look, pli_scores, by="gene") %>% filter(!is.na(lof.pLI))
  missing = sum(!(input.genes %in% matches$gene)) #only 350/375 have scores...
  
  dist <- matches %>% mutate("highly_constrained" = case_when(
    lof.pLI > 0.9 ~ "constrained",
    lof.pLI < 0.1 ~ "not constrained",
    lof.pLI > 0.1 & lof.pLI < 0.9 ~ "indeterminate"
  )
  ) %>% select(gene, gene_id, lof.pLI, highly_constrained,lof.oe,lof.oe_ci.upper,lof.oe_ci.lower)
  
  list("dist" = table(dist$highly_constrained), "missing" = missing, "df"=dist)
}


getPLIDat<- function()
{
  #Do we loose genes b/c of filters?
  pli_scores <- fread("/data/abattle4/aomdahl1/reference_data/gnomad/gnomad.v4.1.constraint_metrics.tsv") %>% filter(!is.na(gene)) %>%
    filter(canonical == TRUE, mane_select == TRUE) %>% group_by(gene) %>% slice_head(n=1)
  
  diffs <- pli_scores %>% group_by(gene) %>% summarize("count" = n(), "diff_lof"= length(unique(lof.pLI)), "range_lof"=max(lof.pLI) - min(lof.pLI))
  stopifnot(!any(diffs$diff_lof != 1))
  #issue.genes <- diffs %>% filter(diff_lof != 1)
  #View(issue.genes %>% filter(range_lof > 0.1))
  #dead.genes <- (issue.genes %>% filter(range_lof > 0.01))$gene
  #View(pli_scores %>% filter(gene %in% dead.genes))
  #If we filter to mane select
  pli_scores
}
