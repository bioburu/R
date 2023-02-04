cat("Please note that any findings in Seurat can be validated with edgeR")
library(dplyr)
library(Seurat)
library(patchwork)
cat("load in 10X data human reference genome version 19 files")
pbmc.data <- Read10X(data.dir = "/home/amplified_prog/RstudioProjects/Seurat/filtered_gene_bc_matrices/hg19")
c("create seurat object with raw counts")
pbmc <- CreateSeuratObject(counts= pbmc.data, project= "pbmc3k", min.cells =3, min.feature=200)
pbmc
cat("scan for counts of molecules in first 30 cells")
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
cat("check memory size of data")
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size/sparse.size
#------------------------------------------------------------------
cat("Remove mitochrondial dna contamination")
cat("Check metadata columns")
head(pbmc)
cat("nCount_RNA= total number of molecules per cell")
cat("nFeature_RNA= total number of genes per cell")
cat("Add new metadata column representing percentage of mito contamination")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc)
#---------------------QC cells--------------------------------------------------
cat("Display plots of QC parameters")
VlnPlot(pbmc, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3)
cat("Produce scatter plots by QC features")
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 =  "nFeature_RNA")
plot2
plot1 + plot2
cat("Data clean up by removing cells with feature counts >2500 or <200")
cat("Data clean up by removing cells with >5% mito counts")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt <5)
head(pbmc)
cat("Visual confirmations of clean up")
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 =  "nFeature_RNA")
plot2
plot1 + plot2
cat("Data is now clean for use")
cat("Normalize counts by dividing gene counts by total gene count and multiplying by 10,000 to scale (similar to cpm in edgeR)")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
head(pbmc[["RNA"]]@data)
cat("All genes are now normalized against each other")
cat("Perform differential expression analysis with genes")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
cat("Get top 10 differentially expressed genes")
top10 <- head(VariableFeatures(pbmc), 10)
top10
dim(pbmc)
plot1 <- VariableFeaturePlot(pbmc)
plot1
#------plot 2 does not seem to be working----------
plot2 <- LabelPoints(plot =plot1, points =top10, repel = TRUE)
plot2
plot1+ plot2
#--------------------------------------------------
cat("Scale all genes to 0 mean with variance of 1")
all.genes <- rownames(pbmc)
all.genes
pbmc <- ScaleData(pbmc, features = all.genes)
head(pbmc[["RNA"]]@scale.data)
cat("All counts are now scaled")
cat("Perform PCA and add to pbmc file")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object=pbmc))
cat("Confirm PCA")
pbmc
cat("Show top 10 DF genes per principle component")
print(pbmc[["pca"]],dims = 1:10,nfeatures = 10)
cat("Visually inspect PCA1 and PCA2")
VizDimLoadings(pbmc,dim=1:2,reduction="pca")
cat("Visually inspect PCA3")
VizDimLoadings(pbmc,dim=3,reduction="pca")
cat("PCA1 vs PCA2")
DimPlot(pbmc, reduction="pca", dims = 1:2)
cat("PCA2 vs PCA3")
DimPlot(pbmc,dim=3:2,reduction="pca")
cat("Inspect first 15 PCA heatmaps")
DimHeatmap(pbmc,dims=1:15,cells=500,balanced = TRUE)
cat("Create null distribution feature scores by randomly permuting a subset of data for rerun of PCA")
pbmc <- JackStraw(pbmc,num.replicate = 100)
cat("Calculate significance scores against null distribution")
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc)
cat("Generate elbow plot")
ElbowPlot(pbmc)
cat("Determine appropriate number of clusters using elbow plot")
cat("K nearest neighbor algorithm")
pbmc <- FindNeighbors(pbmc, dims=1:10)
cat("Generate clusters from KNN results")
pbmc <- FindClusters(pbmc, resolution = 0.5)
cat("Look at first 15 cells and cluster IDs")
head(Idents(pbmc),15)
cat("add UMAP data")
pbmc <- RunUMAP(pbmc,dims=1:10)
DimPlot(pbmc,reduction = "umap", label = TRUE)
pbmc
cat("add TSNE data")
pbmc <- RunTSNE(pbmc,dims=1:10)
DimPlot(pbmc,reduction = "tsne", label = TRUE)
cat("Look at statistics for top 50 genes in cluster 2")
cluster2.markers <- FindMarkers(pbmc,ident.1 = 2, min.pct=0.25)
head(cluster2.markers,n=10)
cat("Search for differences in cluster 2 compared to cluster 0 and 3 only")
cluster2.0_3.markers <- FindMarkers(pbmc,ident.1 = 2,ident.2 = c(0,3),min.pct = 0.25)
head(cluster2.0_3.markers,n=10)
#--------------------?????
pbmc.markers <- FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
pbmc.markers %>% 
  group_by(cluster) %>%
  slice_max(n=2,order_by=avg_log2FC)
#---------------------
cat("Use violin plots to visualize expression patterns")
VlnPlot(pbmc,features = c("CCR7","GZMB","IL7R","IL32","CD3D","S100A9"))
cat("")
cat("Scan features with TSNE")
FeaturePlot(pbmc,features=c("CCR7","GZMB","IL7R","IL32","CD3D","S100A9"),reduction= "tsne")
cat("Scan features with UMAP")
FeaturePlot(pbmc, features =c("CCR7","GZMB","IL7R","IL32","CD3D","S100A9"),reduction = "umap" )
cat("Scan features with PCA")
FeaturePlot(pbmc, features =c("CCR7","GZMB","IL7R","IL32","CD3D","S100A9"),reduction = "pca" )
cat("Add labels to clusters")
new.cluster.ids <- c("Naive CD4 T","CD14+ Mono","Memory CD4T","B","CD8 T","FCGR3A+ Mono",
                     "NK","DC","Platelet")
levels(pbmc)
names(new.cluster.ids) <- levels(pbmc)
new.cluster.ids
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc)
