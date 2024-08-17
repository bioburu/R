library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(htmlwidgets)
require('Seurat.utils')
library(plotly)
#library(caTools)
#library(car)
#library(caret)
#library(InformationValue)
#library(pROC)
#library(ROCR)
fragment<-file.path('/home/deviancedev01/work_stuff/multiome_processing/atac_fragments.tsv.gz')
#---extract barcodes
total_counts <- CountFragments(fragment)
head(total_counts)
summary(total_counts$frequency_count)
cutoff <- 1000 
barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB
head(barcodes)
#-----create fragment object
frags <- CreateFragmentObject(path = fragment, cells = barcodes)
frags
summary(frags@cells)
peaks <- CallPeaks(frags,macs2.path = '/home/deviancedev01/anaconda3/bin/macs3')
peaks
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
head(counts)
#View(data.frame(counts))
head(colnames(counts))
head(row.names(counts))
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = frags,
  min.cells = 10,
  min.features = 200,
  genome = 'hg38')
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'ATAC',
  project = 'scATAC')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotations
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
annotations
Annotation(atac) <- annotations
atac[['ATAC']]
granges(atac)
head(atac)
summary(atac$nCount_ATAC)
VlnPlot(atac,features=c('nCount_ATAC','nFeature_ATAC'),
        pt.size=0.1,
        ncol=2,
        cols = 'skyblue')
atac<-subset(atac,subset=nCount_ATAC>1000)
VlnPlot(atac,features=c('nCount_ATAC','nFeature_ATAC'),
        pt.size=0.1,
        ncol=2,
        cols = 'skyblue')
cat('subset out all TSS scores below 4')
atac <- NucleosomeSignal(object = atac)
atac <- TSSEnrichment(atac)
DensityScatter(atac, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
head(atac)
VlnPlot(object = atac,
        features = c('nCount_ATAC',
                     'nFeature_ATAC',
                     'nucleosome_signal',
                     'nucleosome_percentile',
                     'TSS.enrichment','TSS.percentile'),
        pt.size = 0.1,
        ncol = 3,
        cols = 'skyblue')
atac$high.tss <- ifelse(atac$TSS.enrichment > 4, 'High', 'Low')
table(atac$high.tss)
#TSSPlot(data, group.by = 'high.tss')# + NoLegend()
atac<-subset(atac,subset=TSS.enrichment>4)
VlnPlot(object = atac,
        features = c('nCount_ATAC',
                     'nFeature_ATAC',
                     'nucleosome_signal',
                     'nucleosome_percentile',
                     'TSS.enrichment','TSS.percentile'),
        pt.size = 0.1,
        ncol = 3,
        cols = 'skyblue')
atac <- RunTFIDF(atac)
atac
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac) 
DepthCor(atac)
cat('remove highly correlated dimensions')
atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:30)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims =2:30)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3,resolution = 0.1)
Idents(atac)
atac
DimPlot(atac,label = TRUE)
top500.atac<-FindTopFeatures(atac)[1:500]
row.names(top500.atac)
cat('cleaned ATACseq data is now ok for gene activity analysis')
gene.activities <- GeneActivity(atac)
atac[['correlated_gene_activity']] <- CreateAssayObject(counts = gene.activities)
cat('transfer to gene activity layer')
DefaultAssay(atac)<-'correlated_gene_activity'
atac <- NormalizeData(
  object = atac,
  assay = 'correlated_gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_ATAC))
atac<-FindVariableFeatures(atac)
atac<-ScaleData(atac)
atac
top1000.atac <- head(VariableFeatures(atac), 1000)
top1000.atac
VariableFeatures(atac)
DefaultAssay(atac)<-'ATAC'
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
features_path <- '/home/deviancedev01/work_stuff/multiome_processing/filtered_feature_bc_matrix/features.tsv.gz'
barcodes_path <- '/home/deviancedev01/work_stuff/multiome_processing/filtered_feature_bc_matrix/barcodes.tsv.gz'
matrix_path <- '/home/deviancedev01/work_stuff/multiome_processing/filtered_feature_bc_matrix/matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
colnames(matrix)
row.names(matrix)
rna <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'rna')
summary(rna@active.ident)
gc()
#rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 11000)# & percent.mt <5)
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
rna <- NormalizeData(rna)
gc()
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
VariableFeatures(rna)
rna
rna_df<-data.frame(rna[['RNA']]$counts)
atac_df<-data.frame(rna[['RNA']]$counts)
top5000.rna <- head(VariableFeatures(rna), 5000)
top5000.rna
plot1 <- VariableFeaturePlot(rna)
plot1
all.genes <- rownames(rna)
all.genes
rna <- ScaleData(rna, features = all.genes)
gc()
dim(rna)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))
rna <- RunUMAP(rna, dims = 1:30)
DimHeatmap(rna, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(rna)
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna,resolution = 0.1)
DimPlot(object = rna, label = TRUE,pt.size = 1)
VlnPlot(rna, features = c('PTPRC','CD3D','CD4','CD8A','IFNG','CD19','CD22','CD86','ITGAX'),cols = c())
FeaturePlot(
  object = rna,
  features = c('PTPRC','CD3D','CD4','CD8A','IFNG','CD19','CD22','CD86','ITGAX'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c('grey','red'))
DimPlot(object = rna, label = TRUE,pt.size = 1)
Idents(rna)
levels(rna@meta.data$seurat_clusters)<-c('cytoxic_t_cells','t_helper_cells','mono/macro','b_cells')
rna@meta.data$seurat_clusters
head(rna)
new.cluster.ids <- c('cytotoxic_t_cells','t_helper_cells','mono/macro','b_cells')
names(new.cluster.ids) <- levels(rna)
rna <- RenameIdents(rna, new.cluster.ids)
DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
head(atac)
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "correlated_gene_activity", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$seurat_clusters,
                                     weight.reduction = atac[["lsi"]], dims = 1:30)
atac <- AddMetaData(atac, metadata = celltype.predictions)
head(atac)
head(Idents(atac))
head(Idents(rna))
rna$seurat_clusters
new.cluster.ids <- c('cytotoxic_t_cells','t_helper_cells','mono/macro','b_cells')
names(new.cluster.ids) <- levels(atac)
names(new.cluster.ids)
new.cluster.ids
levels(atac)
Idents(atac)
atac <- RenameIdents(atac, new.cluster.ids)
head(Idents(atac))
head(Idents(rna))
DimPlot(object = atac, label = TRUE,pt.size = 1)
rna_umap <- DimPlot(rna, label = TRUE,pt.size = 1,label.size = 8)+ ggtitle("gene_transcripts")
atac_features<-FeaturePlot(atac,features = c('PTPRC','CD3D','CD4','ICD8A','IFNG','CD19',
                              'CD22','CD86','ITGAX'),cols = c('grey','red'),pt.size = 2)
atac_umap<-DimPlot(object = atac, label = TRUE,pt.size = 1,label.size = 8)
rna_umap+atac_umap
#-------------------------------------------------------------------
peaks <- FindMarkers(
  object = atac,
  ident.1 = 'mono/macro',
  test.use = 'negbinom',
  latent.vars = 'nCount_ATAC',
  logfc.threshold=1.5,
  min.pct=0.2)
peaks
FeaturePlot(atac,features = c(row.names(peaks)[1:12]))
VlnPlot(
  object = atac,
  features = rownames(peaks)[1:12],
  pt.size = 0.1,
  idents = c())
ClosestFeature(atac,c('chr4-8409026-8410072'))
p1<-CoveragePlot(
  object = atac,
  region = 'chr4-8409026-8410072',
  extend.upstream = 40000,
  extend.downstream = 20000)
p2 <- TilePlot(
  object = atac,
  region = 'chr4-8409026-8410072')
p3<-ExpressionPlot(
  object = rna,
  features = c("ACOX3"))
p1
p2+p3
break 
rna
Idents(rna)
data3d <- rna
data3d <- RunUMAP(data3d,dims = 1:15,n.components = 3L)
head(Embeddings(data3d,reduction = "umap"))
rna@meta.data$seurat_clusters

plot3d1 <- FetchData(data3d, vars = c("umap_1", "umap_2", "umap_3", "seurat_clusters"))
plot3d1$label <- paste(plot3d1$seurat_clusters)
fig <- plot_ly(data = plot3d1, 
               x = ~umap_1, y = ~umap_2, z = ~umap_3, 
               color = ~seurat_clusters, 
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
               marker = list(size = 2, width=2),
               text=~label,
               hoverinfo="text")%>%layout(title='scRNAseq_multiome')
fig
#setwd('/home/em_b/work_stuff/snRNA+ATAC_huPBMC')
#saveWidget(ggplotly(fig), file = "multiome_rna.html")

head(data3d@meta.data)
goi <- "CD14"
plotting.data <- FetchData(data3d, vars = c("umap_1", "umap_2", "umap_3","seurat_clusters","Expression"=goi), layer = 'data')
#Cutoff <- quantile(plotting.data[,goi], probs = .95)
#plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")
head(plotting.data)
break
fig2<-plot_ly(data = plotting.data,
        name = goi,
        x = ~umap_1, y = ~umap_2, z = ~umap_3, 
        color = plotting.data$CD14, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 3), 
        text=~seurat_clusters,
        hoverinfo="text"
) %>%layout(title=goi)
#saveWidget(ggplotly(fig2), file = "multiome_rna_cd14.html")
#-------------------------------------------------------------
break 
atac
Idents(atac)
data3d <- atac

head(data3d@meta.data)
data3d <- RunUMAP(object = data3d, reduction = 'lsi',
                  dims = 2:30,
                  n.components = 3L)
levels(atac)
atac@meta.data$predicted.id
#data3d <- ScaleData(data3d)#, features = all.genes)
#data3d <- RunPCA(data3d, features = VariableFeatures(object = top1000.atac))
#data3d
#data3d <- RunUMAP(data3d,dims = 1:15,n.components = 3L)
head(Embeddings(data3d,reduction = "umap"))
plot3d1 <- FetchData(data3d, vars = c("umap_1", "umap_2", "umap_3","TSS.enrichment","predicted.id"))
plot3d1$label<-paste(plot3d1$predicted.id)
fig3 <- plot_ly(data = plot3d1, 
               x = ~umap_1, y = ~umap_2, z = ~umap_3, 
               color = ~predicted.id, 
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
               marker = list(size = 3, width=2),
               text=~label,
               hoverinfo="text")%>%layout(title='scATACseq_multiome')
fig3
#saveWidget(ggplotly(fig3), file = "multiome_atac.html")
fig4 <- plot_ly(data = plot3d1, 
                x = ~umap_1, y = ~umap_2, z = ~umap_3, 
                color = ~TSS.enrichment, 
                colors = c('skyblue','red'),
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 3, width=2),
                text=~label,
                hoverinfo="text")%>%layout(title='scATACseq_multiome')
fig4
#saveWidget(ggplotly(fig4), file = "multiome_atac_tss.html")
break 
#--------------------------------s
saveRDS(rna,file = 'rna.rds')
saveRDS(atac,file='atac.rds')
rna<-readRDS('/home/em_b/Desktop/snRNA+ATAC_huPBMC/rna.rds')
atac<-readRDS('/home/em_b/Desktop/snRNA+ATAC_huPBMC/atac.rds')
