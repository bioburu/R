library(Seurat)
library(readr)
library(Matrix)
library(dplyr)
library(patchwork)
#------------------------------Begin user entries---------------------------------------------------------
#------------------------------Begin user entries---------------------------------------------------------
#------------------------------Begin user entries---------------------------------------------------------
#------------------------------Begin user entries---------------------------------------------------------
#------------------------------Begin user entries---------------------------------------------------------
#Set working directory for experiment
setwd("/home/amp_prog/rstudio/SpeedingCARs")
#Set path to features file
features_path <- "/home/amp_prog/rstudio/SpeedingCARs/D3_hT/features.tsv.gz"
#Set path to barcode file
barcode_path <- "/home/amp_prog/rstudio/SpeedingCARs/D3_hT/barcodes.tsv.gz"
#Set path to matrix file
matrix_path <- "/home/amp_prog/rstudio/SpeedingCARs/D3_hT/matrix.mtx.gz"
#Set minimum number of cells per gene
min_cells <- 50
#Set minimum number of features per cell
min_features <- 200
#Set number of top genes to compare
top_genes <- 2000
#-------------------End User entries and Run with ctrl+shift+enter--------------------------------
#-------------------End User entries and Run with ctrl+shift+enter--------------------------------
#-------------------End User entries and Run with ctrl+shift+enter--------------------------------
#-------------------End User entries and Run with ctrl+shift+enter--------------------------------
#-------------------End User entries and Run with ctrl+shift+enter--------------------------------
cat("Identification of new markers associated with long-term efficacy in patients treated with CAR T cells is a current medical need, particularly in diseases such as multiple myeloma. In this study, we address the impact of CAR density on the functionality of BCMA CAR T cells. Functional and transcriptional studies demonstrate that CAR T cells with high expression of the CAR construct show an increased tonic signaling with up-regulation of exhaustion markers and increased in vitro cytotoxicity but a decrease in in vivo BM infiltration. Characterization of gene regulatory networks using scRNA-seq identified regulons associated to activation and exhaustion up-regulated in CARHigh T cells, providing mechanistic insights behind differential functionality of these cells. Last, we demonstrate that patients treated with CAR T cell products enriched in CARHigh T cells show a significantly worse clinical response in several hematological malignancies. In summary, our work demonstrates that CAR density plays an important role in CAR T activity with notable impact on clinical response.")
cat("Displaying features file")
read_tsv(features_path)
cat("Displaying barcodes file")
read_tsv(barcode_path)
cat("Displaying expression matrix")
readMM(matrix_path)
cat("Creating read matrix with labels")
matrix <- ReadMtx(mtx = matrix_path, features = features_path, cells = barcode_path)
print(matrix)
dim(matrix)
cat("Create a S4 object")
object <- CreateSeuratObject(counts=matrix, min.cells = min_cells, min.features = min_features)
rm(matrix)
gc()
dim(object)
head(object)
cat("Scan for mitochondrial contamination")
object[["percent_mito"]] <- PercentageFeatureSet(object, pattern = "^MT-")
head(object)
cat("Displaying QC parameters")
VlnPlot(object, features = c("nCount_RNA","nFeature_RNA","percent_mito"),ncol = 3)
FeatureScatter(object, feature1 = "nCount_RNA",feature2 = "percent_mito")
FeatureScatter(object, feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
cat("Removing feature counts >5500 or <200")
cat("Removing >5% mitochrondial contamination")
car_t <- subset(object, subset= nFeature_RNA >200 & nFeature_RNA <6000 & percent_mito <5)
head(car_t)
rm(object)
gc()
cat("Data is now clean for use")
cat("Normalize ie reads/total_reads*10,000")
car_T <- NormalizeData(car_t, normalization.method = "LogNormalize", scale.factor = 10000)
head(car_T[["RNA"]]@data)
rm(car_t)
gc()
cat("All genes are now normalized")
cat("Isolate most variable genes")
Car_T <- FindVariableFeatures(car_T, selection.method = "vst", nfeatures = top_genes)
head(Car_T)
rm(car_T)
gc()
cat("Display top 200 genes")
head(VariableFeatures(Car_T), 200) 
cat("Display Variance plot")
VariableFeaturePlot(Car_T)
cat("Scale all reads to 0")
all.genes <- rownames(Car_T) 
print(all.genes)
gc()
CAr_T <- ScaleData(Car_T, features = all.genes)
head(CAr_T)
rm(Car_T)
gc()
cat("Perform PCA")
CAr_T <- RunPCA(CAr_T, features = VariableFeatures(object=CAr_T))
print(CAr_T)
cat("Top 10 components")
print(CAr_T[["pca"]], dim=1:10,nFeatures=10)
cat("PCA1 PCA2")
VizDimLoadings(CAr_T,dims = 1:2,reduction="pca")
cat("PCA2 PCA3")
VizDimLoadings(CAr_T,dims = 2:3,reduction = "pca")
cat("PCA1vPCA2")
DimPlot(CAr_T,reduction = "pca",dims = 1:2)
cat("PCA2vPCA3")
DimPlot(CAr_T,reduction = "pca",dims = 2:3)
cat("PCA1-15")
DimHeatmap(CAr_T,dims = 1:15,cells=500,balanced = TRUE)
cat("PCA16-30")
DimHeatmap(CAr_T,dims = 16:30,cells = 500,balanced = TRUE)
cat("PCA31-45")
DimHeatmap(CAr_T,dims = 31:45,cells = 500,balanced = TRUE)
cat("Create null distribution")
gc()
rm(all.genes,barcode_path,features_path,min_cells,min_features,top_genes,matrix_path)
gc()
CAr_T <- JackStraw(CAr_T,num.replicate = 100)
gc()
CAR_T <- ScoreJackStraw(CAr_T,reduction= "pca", dims = 1:20)
JackStrawPlot(CAR_T)
cat("Significant scores created")
rm(CAr_T)
gc()
cat("Generating elbow")
ElbowPlot(CAR_T,ndims = 50)
cat("K nearest neighbor algorithm")
CAR_T <- FindNeighbors(CAR_T)
CAR_T
cat("Generating clusters")
CAR_T <- FindClusters(CAR_T,resolution = 1.0)
head(Idents(CAR_T),15)
cat("Enter in number of dimensions to use for UMAP in console o/")
num=as.integer(readline(prompt = "Enter number of clusters for dimensional reductions: "))
CAR_T <- RunUMAP(CAR_T, dims = 1:num)
DimPlot(CAR_T,reduction = "umap",label = TRUE)
gc()
cat("Creating t-SNE")
CAR_T <- RunTSNE(CAR_T,dims = 1:num)
DimPlot(CAR_T,reduction = "tsne",label = TRUE)
cat("Statistics for cluster 0")
head(FindMarkers(CAR_T,ident.1 = 0,min.pct = 0.25),n=50)
cat("Statistics for cluster 1")
head(FindMarkers(CAR_T,ident.1 = 1,min.pct = 0.25),n=50)
cat("Statistics for cluster 2")
head(FindMarkers(CAR_T,ident.1 = 2,min.pct = 0.25),n=50)
gc()
cat("Statistics for cluster 3")
head(FindMarkers(CAR_T,ident.1 = 3,min.pct = 0.25),n=50)
cat("Statistics for cluster 4")
head(FindMarkers(CAR_T,ident.1 = 4,min.pct = 0.25),n=50)
gc()
cat("Statistics for cluster 5")
head(FindMarkers(CAR_T,ident.1 = 5,min.pct = 0.25),n=50)
gc()
#markers <- FindAllMarkers(CAR_T,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
#markers %>% 
#  group_by(cluster) %>%
#  slice_max(n=2,order_by=avg_log2FC)
cat("Select and enter genes of interest")
break
#-------------Begin user entries for genes----------------------------------------
#-------------Begin user entries for genes----------------------------------------
#-------------Begin user entries for genes----------------------------------------
#-------------Begin user entries for genes----------------------------------------
#-------------Begin user entries for genes----------------------------------------
#-------------Begin user entries for genes----------------------------------------
#Enter genes for violin plots
VlnPlot(CAR_T,features = c("","","","","",""))
#Enter genes for tSNE        
FeaturePlot(CAR_T,features=c("","","","","",""),reduction= "tsne")
#Enter genes for umap
FeaturePlot(CAR_T,features=c("","","","","",""),reduction= "umap")
#Enter genes for pca
FeaturePlot(CAR_T,features=c("","","","","",""),reduction= "pca")
#All done \o/
