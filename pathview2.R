#--------------
library(pathview)
foldchanges = df$log2FoldChange
names(foldchanges) = df$entrezgene_id
head(foldchanges)
setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/kegg')
pv_bulk_tissues <- pathview(gene.data = foldchanges, pathway.id = "04060",
                   species = "mmu", out.suffix = "mmu_test")
 
