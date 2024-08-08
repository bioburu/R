library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(plotly)
library(htmlwidgets)
setwd('/home/em_b/work_stuff/differentiation_bulkrna')
matrix<-read.csv('rnaseq_bulk_0hr_48hr_48t3.csv')
row.names(matrix)<-make.names(matrix$external_gene_name,
                              unique = TRUE)
matrix2<-matrix[,-c(1:2)]
matrix2<-round(matrix2)
data <- CreateSeuratObject(counts=matrix2,project = 'MB_GNP.T3')
gc()
table(data@active.ident)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(data))
top1000 <- head(VariableFeatures(data), 1000)
head(top1000)
gc()
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
length(all.genes)
gc()
data <- ScaleData(data, features = all.genes)
dim(data)
data <- RunPCA(data,npcs = 17, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
table(data@meta.data$orig.ident)
#data <-JoinLayers(data)
data3d<-data.frame(data@reductions$pca@cell.embeddings)
str(data3d)
#plot3d <- FetchData(data3d, vars = c("PC_1", "PC_2", "PC_3", row.names()))
data3d$label <- paste(rownames(data3d))
data3d$label <-gsub("\\_.*","",data3d$label)
fig <- plot_ly(data = data3d, 
               x = ~PC_9, y = ~PC_14, z = ~PC_1, 
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
               marker = list(size = 10, width=5),
               text=~label,
               hoverinfo="text")
#----------------------------------------------------
df<-data.frame(t(data@assays$RNA$counts))
df$Neurod1
data3d2<-cbind(df$Neurod1,data3d)
Neurod1<-data3d2$`df$Neurod1`
fig2 <- plot_ly(data = data3d, 
               x = ~PC_9, y = ~PC_14, z = ~PC_1, 
               color = ~Neurod1, 
               colors = c(),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 10, width=5),
               text=data3d$label,
               hoverinfo="text")%>%layout(title='Neurod1')
#----------------------------------------------------
df<-data.frame(t(data@assays$RNA$counts))
df$Gli2
data3d2<-cbind(df$Gli2,data3d)
Gli2<-data3d2$`df$Gli2`
fig3 <- plot_ly(data = data3d, 
                x = ~PC_9, y = ~PC_14, z = ~PC_1, 
                color = ~Gli2, 
                colors = c(),
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 10, width=5),
                text=data3d$label,
                hoverinfo="text")%>%layout(title='Gli2')
#----------------------------------------------------
df<-data.frame(t(data@assays$RNA$counts))
df$Hr
data3d2<-cbind(df$Hr,data3d)
Hr<-data3d2$`df$Hr`
fig4 <- plot_ly(data = data3d, 
                x = ~PC_9, y = ~PC_14, z = ~PC_1, 
                color = ~Hr, 
                colors = c(),
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 10, width=5),
                text=data3d$label,
                hoverinfo="text")%>%layout(title='Hr')
#----------------------------------------------------
df<-data.frame(t(data@assays$RNA$counts))
df$Cntn2
data3d2<-cbind(df$Cntn2,data3d)
Cntn2<-data3d2$`df$Cntn2`
fig5 <- plot_ly(data = data3d, 
                x = ~PC_9, y = ~PC_14, z = ~PC_1, 
                color = ~Cntn2, 
                colors = c(),
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 10, width=5),
                text=data3d$label,
                hoverinfo="text")%>%layout(title='Cntn2')
#----------------------------------------------------
df<-data.frame(t(data@assays$RNA$counts))
df$Rora
data3d2<-cbind(df$Rora,data3d)
Rora<-data3d2$`df$Rora`
fig6 <- plot_ly(data = data3d, 
                x = ~PC_9, y = ~PC_14, z = ~PC_1, 
                color = ~Rora, 
                colors = c(),
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 10, width=5),
                text=data3d$label,
                hoverinfo="text")%>%layout(title='Rora')
#----------------------------------------------------
df<-data.frame(t(data@assays$RNA$counts))
df$Gli1
data3d2<-cbind(df$Gli1,data3d)
Gli1<-data3d2$`df$Gli1`
fig7 <- plot_ly(data = data3d, 
                x = ~PC_9, y = ~PC_14, z = ~PC_1, 
                color = ~Gli1, 
                colors = c(),
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 10, width=5),
                text=data3d$label,
                hoverinfo="text")%>%layout(title='Gli1')



break
fig
fig2
fig3
fig4
fig5
fig6
fig7
p <- plotly_build(fig1)
p$x$data[[3]]$marker$symbol <- 'diamond'
p$x$data[[6]]$marker$symbol <- 'diamond'
p
saveWidget(ggplotly(fig7), file = "Gli1.html")
break 
#-----correlations
FeatureScatter(data,
               feature1 = "Neurod1",
               feature2 = 'Ezh2',
               cols = c(),
               pt.size = 5,
               shuffle = TRUE,
               seed = 1,
               jitter = TRUE)
VlnPlot(data,
        features = c(),
        pt.size = 1,
        cols = c('grey','red',
                 'grey','red'),
        layer='counts')

DimPlot(data,
        cols = c('grey','black','red',
                 'grey','black','red'),
        pt.size = 7,
        dims = c(14,9))
#-----------PCA features
TopFeatures(data[['pca']],
            dim = 9,
            balanced = TRUE)
#------ridgeplot
RidgePlot(data,
          features = c(),
          cols = c())
#------------------------------------------
x<-FindMarkers(data, ident.1 = 'GBM', ident.2 = 'NSC', 
               features = c(top1000),logfc.threshold=2,
               only.pos = FALSE,test.use = 'DESeq2',min.pct = 0.5)
VlnPlot(data,features = c(rownames(x)[1:12]),cols = c())
x<-subset(x,p_val_adj< 0.01)
downreg<-subset(x,avg_log2FC< -1.5)
upreg<-subset(x,avg_log2FC> 1.5)
df<-rbind(upreg,downreg)
df<-cbind(row.names(df),df)
break
write.csv(x,file = 'GSE154958_gbm_dura.csv')

