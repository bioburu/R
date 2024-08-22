library(Seurat)
library(ShinyCell)
setwd('/home/deviancedev01/Desktop')
data<-readRDS('/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/rna_p7.rds')
scConf<-createConfig(data)
makeShinyApp(data,
             scConf,
             gene.mapping = TRUE,
             shiny.title = 'MB_tumors_ptch1.mth1_p7_npcs')
