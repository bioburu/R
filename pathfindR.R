BiocManager::install("org.Hs.eg.db")
#-------------
sudo apt-get install libmagick++-dev
#------------
install.packages("magick")
BiocManager::install("KEGGgraph")
install.packages('pathfindR')
library(pathfindR)
#----obtain mouse specific gene sets
#gsets_list<-get_gene_sets_list(source = 'mmu_KEGG')
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/alignments/finals')
matrix<-read.csv('DESeq2_out.csv')
matrix<-matrix[order(matrix$p_val_adj), ]
matrix<-matrix[,-c(1,3,5:7)]
#------clean gene set
#genes<-input_processing(input = matrix,p_val_threshold = 0.01,pin_name_path = 'KEGG')
#------get gene sets
gene_list <- fetch_gene_set(
  gene_sets = "mmu_KEGG",
  min_gset_size = 10,
  max_gset_size = 300)
#------------------------------------------------------------------------
mmu_kegg_gsets <- gene_list[[1]]
head(mmu_kegg_gsets)
#------------------------------------------------------------
mmu_kegg_descriptions <- gene_list[[2]]
head(mmu_kegg_descriptions)
#---
break
matrix<-matrix[,-2]
test<-run_pathfindR(
  input = matrix,
  convert2alias = FALSE,
  custom_genes = mmu_KEGG_gsets,
  custom_descriptions = mmu_KEGG_descriptions,
  pin_name_path = 'mmu_STRING'
)
cluster_enriched_terms(test)
visualize_terms(test)
