BiocManager::install("org.Hs.eg.db")
#-------------
sudo apt-get install libmagick++-dev
#------------
install.packages("magick")
BiocManager::install("KEGGgraph")
install.packages('pathfindR')
library(pathfindR)
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/alignments/finals/results')
matrix<-read.csv('DESeq2_out.csv')
#matrix<-matrix[order(matrix$p_val_adj), ]
matrix<-matrix[,-c(1,3,5)]
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
  convert2alias = FALSE,
  custom_genes = mmu_KEGG_gsets,
  custom_descriptions = mmu_KEGG_descriptions,
  pin_name_path = 'mmu_STRING'
)
data<-cluster_enriched_terms(test)
term_gene_graph(test)
term_gene_graph(test,num_terms = c(5),use_description = TRUE)
enrichment_chart(test, plot_by_cluster = TRUE)
data<-data[order(data$highest_p), ]
#write.csv(data,file = 'mmu_KEGG_results.csv')
break 
#----for individual cluster plots
visualize_terms(result_df = data,
                hsa_KEGG = FALSE,
                pin_name_path = 'mmu_STRING')


