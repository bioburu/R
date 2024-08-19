library(Signac)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(htmlwidgets)
require('Seurat.utils')
library(plotly)
library(EnsDb.Mmusculus.v79)
library(celldex)
library(SingleR)
setwd('/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis')
#-------------------------------------------------------------------------------
features_path <- '/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/p7_gnp/outs/filtered_feature_bc_matrix/features.tsv.gz'
barcodes_path <- '/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/p7_gnp/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
matrix_path <- '/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/p7_gnp/outs/filtered_feature_bc_matrix/matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
#head(colnames(matrix))
#head(row.names(matrix))
rna <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'rna')
summary(rna@active.ident)
#rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
#VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 11000)# & percent.mt <5)
#VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
#head(VariableFeatures(rna))
#rna
#rna_df<-data.frame(rna[['RNA']]$counts)
#atac_df<-data.frame(atac[['ATAC']]$counts)
#top5000.rna <- head(VariableFeatures(rna), 5000)
#top5000.rna
plot1 <- VariableFeaturePlot(rna)
plot1
all.genes <- rownames(rna)
#all.genes
rna <- ScaleData(rna, features = all.genes)
#dim(rna)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))
rna <- RunUMAP(rna, dims = 1:30)
DimHeatmap(rna, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(rna)
rna <- FindNeighbors(rna, dims = 1:30)
#--adjust resolution to atac umap clusters
rna <- FindClusters(rna,resolution = 0.10)
DimPlot(object = rna, label = TRUE,pt.size = 1)
#--auto annotate cell clusters
surveyReferences()
searchReferences('mouse')
ref <- fetchReference("mouse_rnaseq", "2024-02-26")
ref
#--- convert seurat object into single cell experiment object 
rna.sce <- as.SingleCellExperiment(DietSeurat(rna))
rna.sce
ref.main <- SingleR(test = rna.sce,
                       assay.type.test = 1,
                       ref = ref,
                       labels = ref$label.main)
ref.fine <- SingleR(test = rna.sce,
                    assay.type.test = 1,
                    ref = ref,
                    labels = ref$label.fine)
table(ref.main$pruned.labels)
table(ref.fine$pruned.labels)
rna@meta.data$ref.main <- ref.main$pruned.labels
rna@meta.data$ref.fine <- ref.fine$pruned.labels
head(rna@meta.data)
rna <- SetIdent(rna, value = "ref.fine")
FeaturePlot(
  object = rna,
  features = c('Mycn'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c('grey','red'))
FeaturePlot(
  object = rna,
  features = c('Rbfox3'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c('grey','red'))
FeaturePlot(
  object = rna,
  features = c('Gli2'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c('grey','red'))
table(Idents(rna))
#levels(rna@meta.data$seurat_clusters)<-c('MB_tumors','mature_neurons','neurons1','early_neurons',
#                                         'oligodendrocytes','neurons2','microglia','unknown1','unknown2')
#rna@meta.data$seurat_clusters
#head(rna)
#new.cluster.ids <- c('MB_tumors','mature_neurons','neurons1','early_neurons',
#                     'oligodendrocytes','neurons2','microglia','unknown1','unknown2')
levels(rna)
#names(nrna.scenames(new.cluster.ids) <- levels(rna)
#rna <- RenameIdents(rna, new.cluster.ids)
DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
break 
#--start atac here
fragment<-file.path('/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/p7_gnp/outs/atac_fragments.tsv.gz')
#---extract barcodes
total_counts <- CountFragments(fragment)
#head(total_counts)
#summary(total_counts$frequency_count)
cutoff <- 1000 
barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB
#head(barcodes)
#-----create fragment object
frags <- CreateFragmentObject(path = fragment, cells = barcodes)
#frags
#summary(frags@cells)
peaks <- CallPeaks(frags,macs2.path = '/home/deviancedev01/anaconda3/bin/macs3')
#peaks
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
#head(counts)
#View(data.frame(counts))
#head(colnames(counts))
#head(row.names(counts))
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = frags,
  min.cells = 10,
  min.features = 200,
  genome = 'mm10')
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'ATAC',
  project = 'scATAC')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#annotations
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
annotations
Annotation(atac) <- annotations
atac[['ATAC']]
#granges(atac)
#head(atac)
#summary(atac$nCount_ATAC)
#VlnPlot(atac,features=c('nCount_ATAC','nFeature_ATAC'),
#        pt.size=0.1,
#        ncol=2,
#        cols = 'skyblue')
atac<-subset(atac,subset=nCount_ATAC>1000)
#VlnPlot(atac,features=c('nCount_ATAC','nFeature_ATAC'),
#        pt.size=0.1,
#        ncol=2,
#        cols = 'skyblue')
#cat('subset out all TSS scores below 4')
atac <- NucleosomeSignal(object = atac)
atac <- TSSEnrichment(atac)
DensityScatter(atac, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = NULL)
#head(atac)
#VlnPlot(object = atac,
#        features = c('nCount_ATAC',
#                     'nFeature_ATAC',
#                     'nucleosome_signal',
#                     'nucleosome_percentile',
#                     'TSS.enrichment','TSS.percentile'),
#        pt.size = 0.1,
#        ncol = 3,
#        cols = 'skyblue')
atac$high.tss <- ifelse(atac$TSS.enrichment > 4, 'High', 'Low')
table(atac$high.tss)
#TSSPlot(data, group.by = 'high.tss')# + NoLegend()
atac<-subset(atac,subset=TSS.enrichment>4)
VlnPlot(object = atac,
        features = c('nCount_ATAC',
                     'nFeature_ATAC',
                     'nucleosome_signal',
                     'nucleosome_percentile',
                     'TSS.enrichment',
                     'TSS.percentile'),
        pt.size = 0.1,
        ncol = 3,
        cols = 'skyblue')
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac) 
DepthCor(atac)
#-- remove highly correlated dimensions
atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:30)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims =2:30)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3,resolution = 0.1)
#Idents(atac)
#atac
DimPlot(atac,label = TRUE)
#top500.atac<-FindTopFeatures(atac)[1:500]
#row.names(top500.atac)
#-- cleaned ATACseq data is now ok for gene activity analysis'
gene.activities <- GeneActivity(atac)
#View(data.frame(gene.activities))
atac[['correlated_gene_activity']] <- CreateAssayObject(counts = gene.activities)
head(atac@meta.data)
#-- transfer to gene activity layer'
DefaultAssay(atac)<-'correlated_gene_activity'
atac <- NormalizeData(
  object = atac,
  assay = 'correlated_gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_ATAC))
atac<-FindVariableFeatures(atac)
atac<-ScaleData(atac)
#atac
#top1000.atac <- head(VariableFeatures(atac), 1000)
#head(top1000.atac)
#head(VariableFeatures(atac))
DefaultAssay(atac)<-'ATAC'
#atac
#----
#head(atac)
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "correlated_gene_activity", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$ref.fine,
                                     weight.reduction = atac[["lsi"]], dims = 1:30)
atac <- AddMetaData(atac, metadata = celltype.predictions)
#head(atac)
#head(Idents(atac))
#head(Idents(rna))
#rna$seurat_clusters
#new.cluster.ids <- c('MB_tumors','mature_neurons','neurons1','early_neurons',
#                     'oligodendrocytes','neurons2','microglia','unknown1','unknown2')
Idents(atac)<-atac@meta.data$predicted.id
#Idents(atac)
#names(new.cluster.ids) <- levels(atac)
#names(new.cluster.ids)
#new.cluster.ids
#levels(atac)
#Idents(atac)
#atac <- RenameIdents(atac, new.cluster.ids)
#head(Idents(atac))
#head(Idents(rna))
DimPlot(object = atac, label = TRUE,pt.size = 1)
rna_umap <- DimPlot(rna, label = TRUE,pt.size = 1,label.size = 3)+ ggtitle("gene_transcripts")
atac_features<-FeaturePlot(atac,features = c('Sox2','Olig2','Mycn','Mrc1','Rbfox3',
                                             'Mki67','Gli2'),cols = c('grey','red'),pt.size = 1)
atac_features
atac_umap<-DimPlot(object = atac, label = TRUE,pt.size = 1,label.size = 3)+ ggtitle("chromatin_accessibility")
rna_umap+atac_umap
#-------------------------------------------------------------------
#peaks <- FindMarkers(
#  object = atac,
#  ident.1 = 'Neurons',
#  test.use = 'negbinom',
#  latent.vars = 'nCount_ATAC',
#  logfc.threshold=1,
#  min.pct=0.2)
#peaks
#FeaturePlot(atac,features = c(row.names(peaks)[1:12]))
#VlnPlot(
#  object = atac,
#  features = rownames(peaks)[1:12],
#  pt.size = 0.1,
#  idents = c())
ClosestFeature(atac,c('chr1-118834061-119054405'))
p1<-CoveragePlot(
  object = atac,
  region = 'chr1-118834061-119054405',
  extend.upstream = 40000,
  extend.downstream = 20000)
p2 <- TilePlot(
  object = atac,
  region = 'chr1-118834061-119054405')
#p3<-ExpressionPlot(
#  object = rna,
#  features = c("Gli2"))
p1

p2
#p3
#---- 3d plot rna 
data3d <- rna
data3d <- RunUMAP(data3d,dims = 1:15,n.components = 3L)
#head(Embeddings(data3d,reduction = "umap"))
plot3d1 <- FetchData(data3d, vars = c("umap_1", "umap_2", "umap_3", "ref.fine"))
plot3d1$label <- paste(plot3d1$ref.fine)
fig <- plot_ly(data = plot3d1, 
               x = ~umap_1, y = ~umap_2, z = ~umap_3, 
               color = ~ref.fine, 
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
#setwd('/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/p7_gnp')
#saveWidget(ggplotly(fig), file = "multiome_scrnaseq.html")
goi <- "Gli2"
plotting.data <- FetchData(data3d, vars = c("umap_1", "umap_2", "umap_3","ref.fine","Expression"=goi), layer = 'data')
#Cutoff <- quantile(plotting.data[,goi], probs = .95)
#plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)
#plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")
#head(plotting.data)
fig2<-plot_ly(data = plotting.data,
        name = goi,
        x = ~umap_1, y = ~umap_2, z = ~umap_3, 
        color = plotting.data$Gli2, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 3), 
        text=~ref.fine,
        hoverinfo="text"
) %>%layout(title=goi)
fig2
#saveWidget(ggplotly(fig2), file = "multiome_scrnaseq_gli2.html")
#-- 3d plot atac
data3d <- atac
#head(data3d@meta.data)
data3d <- RunUMAP(object = data3d, reduction = 'lsi',
                  dims = 2:30,
                  n.components = 3L)
levels(atac)
#atac@meta.data$predicted.id
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
#saveWidget(ggplotly(fig3), file = "multiome_scatacseq.html")
#fig4 <- plot_ly(data = plot3d1, 
#                x = ~umap_1, y = ~umap_2, z = ~umap_3, 
#                color = ~TSS.enrichment, 
#                colors = c('grey','red'),
#                type = "scatter3d", 
#                mode = "markers", 
#                marker = list(size = 3, width=2),
#                text=~label,
#                hoverinfo="text")%>%layout(title='scATACseq_multiome')
#fig4
#saveWidget(ggplotly(fig4), file = "multiome_scatacseq_tss.html")
#--------------------------------s
#saveRDS(rna,file = 'rna_p7.rds')
#saveRDS(atac,file='atac_p7.rds')
