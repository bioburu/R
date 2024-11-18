library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(plotly)
library(htmlwidgets)
library(DESeq2)
library(volcano3D)
library(org.Mm.eg.db)
data <- CreateSeuratObject(counts=matrix,project = 'carM')
table(data@active.ident)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data,npcs = 8, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:8, cells = 500, balanced = T)
ElbowPlot(data)
table(data@meta.data$orig.ident)
#---extract PCAs for plotting
data3d<-data.frame(data@reductions$pca@cell.embeddings)
head(data3d)
data3d$label <- paste(rownames(data3d))
data3d$label <-gsub("\\_.*","",data3d$label)
fig <- plot_ly(data = data3d, 
               x = ~PC_1, y = ~PC_2, z = ~PC_3, 
               color = ~label, 
               colors = c("lightseagreen",
                          "gray50",
                          "darkgreen",
                          "red4",
                          "red",
                          "turquoise4",
                          "black",
                          "yellow4",
                          "royalblue1",
                          "lightcyan3",
                          "peachpuff3",
                          "khaki3",
                          "gray20",
                          "orange2",
                          "royalblue4",
                          "yellow3",
                          "gray80",
                          "darkorchid1",
                          "lawngreen",
                          "plum2",
                          "darkmagenta"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 15, width=5),
               text=~label,
               hoverinfo="text")%>%layout(title='GSE201174:RNA-seq data for CAR_M and WT_M cells with and without CHLA-neuroblastoma co-culture')
fig
#----------3d plot by gene expressions 
genes<-data.frame(t(data@assays$RNA$counts))
list<-c('IL1B','IL10','CXCL8')
genes<-genes[,c(list)]
data3d<-cbind(data3d,genes)
fig2 <- plot_ly(data = data3d, 
                x = ~IL1B, y = ~IL10, z = ~CXCL8, 
                color = ~label, 
                colors = c(),
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 10, width=5),
                text=data3d$label,
                hoverinfo="text")%>%layout(title='GSE201174:RNA-seq data for CAR_M and WT_M cells with and without CHLA-neuroblastoma co-culture')
fig2
VlnPlot(data,
          features = c('IL1B','IL10','CXCL8'),
          cols = c())
