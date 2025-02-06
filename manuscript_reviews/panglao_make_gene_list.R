#Make panglao file for permutation tests:
pacman::p_load(magrittr, dplyr, ggplot2, data.table)
all.markers <- fread("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
#extract data for permutation tests later:
keep <- c("Platelets","Hematopoietic stem cells","Megakaryocytes")
out <- all.markers %>% filter(`cell type` %in% keep) %>% select(`cell type`,`official gene symbol`) %>% group_by(`cell type`) %>%
  summarize("genes"=paste(`official gene symbol`,collapse = ",")) %>% mutate("cell_type"=gsub(x=`cell type`,pattern = " ",replacement = "_")) %>%
  select(cell_type,genes)
write.table(out, file="/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/panglao_genes.txt", quote = FALSE,row.names = FALSE, col.names = FALSE)
