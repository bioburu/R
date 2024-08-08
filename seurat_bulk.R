library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(plotly)
library(htmlwidgets)
library(DESeq2)
library(volcano3D)
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
#------3d volcancos
gnps<-matrix2[,c(1:9)]
head(gnps)
names<-colnames(gnps)
names
condition<-c('GNP.0hr','GNP.0hr','GNP.0hr',
             'GNP.48hr','GNP.48hr','GNP.48hr',
             'GNP.48hr.T3','GNP.48hr.T3','GNP.48hr.T3')
type<-c('paired','paired','paired',
        'paired','paired','paired',
        'paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
coldata
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
colnames(gnps)
row.names(coldata)

deseq<-DESeqDataSetFromMatrix(countData = gnps,
                              colData = coldata,
                              design = ~ condition)
deseq
DE<-DESeq(deseq)
DE
LRT<-DESeq(deseq,
           test = c('LRT'),
           reduced = ~ 1,
           parallel = TRUE)
res<-deseq_polar(DE,LRT,'condition')
str(res)
vol3d<-volcano3D(
  res,
  type = 2,
  label_rows = c('Neurod1','Gli2','Klf9'),
  label_size = 14,
  arrow_length = 120,
  colour_code_labels = TRUE,
  label_colour = "black",
  grid_colour = "grey80",
  grid_width = 2,
  grid_options = NULL,
  axis_colour = "red",
  axis_width = 2,
  marker_size = 5,
  marker_outline_width = 1,
  marker_outline_colour = "grey",
  z_axis_title_offset = 1.2,
  z_axis_title_size = 12,
  z_axis_angle = 0.5,
  radial_axis_title_size = 20,
  radial_axis_title_offset = 1.2,
  xy_aspectratio = 3,
  z_aspectratio = 1,
  camera_eye = list(x = 0.9, y = 0.9, z = 0.9))
#-----------------
mbs<-matrix2[,c(10:18)]
head(mbs)
names<-colnames(mbs)
names
condition<-c('MB.0hr','MB.0hr','MB.0hr',
             'MB.48hr','MB.48hr','MB.48hr',
             'MB.48hr.T3','MB.48hr.T3','MB.48hr.T3')
type<-c('paired','paired','paired',
        'paired','paired','paired',
        'paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
coldata
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
colnames(mbs)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = mbs,
                              colData = coldata,
                              design = ~ condition)
deseq
DE<-DESeq(deseq)
DE
LRT<-DESeq(deseq,
           test = c('LRT'),
           reduced = ~ 1,
           parallel = TRUE)
res<-deseq_polar(DE,LRT,'condition')
str(res)
vol3d<-volcano3D(
  res,
  type = 2,
  label_rows = c('Neurod1','Gli2','Klf9'),
  label_size = 14,
  arrow_length = 120,
  colour_code_labels = TRUE,
  label_colour = "black",
  grid_colour = "grey80",
  grid_width = 2,
  grid_options = NULL,
  axis_colour = "red",
  axis_width = 2,
  marker_size = 5,
  marker_outline_width = 1,
  marker_outline_colour = "grey",
  z_axis_title_offset = 1.2,
  z_axis_title_size = 12,
  z_axis_angle = 0.5,
  radial_axis_title_size = 20,
  radial_axis_title_offset = 1.2,
  xy_aspectratio = 3,
  z_aspectratio = 1,
  camera_eye = list(x = 0.9, y = 0.9, z = 0.9))
vol3d
break 
#---extract PCAs for plotting
data3d<-data.frame(data@reductions$pca@cell.embeddings)
str(data3d)
head(data3d)
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
fig
#-----------PCA features
TopFeatures(data[['pca']],
            dim = 1,
            balanced = TRUE)
#----------3d plot by gene expressions 
genes<-data.frame(t(data@assays$RNA$counts))
list<-c('Neurod1','Gli2','Hr')
genes<-genes[,c(list)]
data3d<-cbind(data3d,genes)
fig2 <- plot_ly(data = data3d, 
               x = ~Neurod1, y = ~Gli2, z = ~Hr, 
               color = ~label, 
               colors = c(),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 10, width=5),
               text=data3d$label,
               hoverinfo="text")%>%layout(title='T3_treatments_vs_controls')
fig2
break 
saveWidget(ggplotly(fig2), file = "neurod1.gli1.hr.html")
fig
fig2
fig3
fig4
fig5
fig6
fig7
vol3d
break 
#---change shape of plot if needed
p <- plotly_build(fig1)
p$x$data[[3]]$marker$symbol <- 'diamond'
p$x$data[[6]]$marker$symbol <- 'diamond'
p
saveWidget(ggplotly(vol3d), file = "mb.3dvol.html")
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
