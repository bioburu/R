BiocManager::install("org.Hs.eg.db")
#-------------
sudo apt-get install libmagick++-dev
#------------
install.packages("magick")
BiocManager::install("KEGGgraph")
install.packages('pathfindR')
library(pathfindR)
library(dplyr)
setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results')
x<-read.csv('DESeq2_MB.bulk.csv',row.names = 1)
x<-x[order(x$Gene),]
x<-x%>%group_by(Gene)%>%arrange(padj)%>%slice_min(padj)
x<-x[-c(310:358),]
#matrix<-matrix[order(matrix$p_val_adj), ]
matrix<-data.frame(x[,-c(1,3,5:7)])
summary(matrix)
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
#matrix<-matrix[,-2]
test<-run_pathfindR(
  input = matrix,
  gene_sets = 'mmu_KEGG',
  convert2alias = FALSE,
  custom_genes = mmu_KEGG_gsets,
  custom_descriptions = mmu_KEGG_descriptions,
  pin_name_path = 'mmu_STRING'
)
str(matrix)
data<-cluster_enriched_terms(test)
term_gene_graph(test)
term_gene_graph(test,num_terms = c(5),use_description = TRUE)
enrichment_chart(test, plot_by_cluster = TRUE,top_terms = 20,num_bubbles = 8)
data<-data[order(data$highest_p), ]
#write.csv(data,file = 'bulk.tiss_mmu_KEGG_results.csv')
break 
#----for individual cluster plots
visualize_terms(result_df = data,
                hsa_KEGG = FALSE,
                pin_name_path = 'mmu_STRING')
devtools::install_github("egeulgen/pathfindR")
input_processed <- input_processing(matrix)
visualize_terms(
  result_df = test,
  input_processed = input_processed,
  hsa_KEGG = TRUE
)
#----------------------
#fishers exact test can be performed to confirm cluster enrichements using a contingency table 
#https://carpentries-incubator.github.io/bioc-rnaseq/06-gene-set-analysis.html
