library(Seurat)
library(readr)
library(Matrix)
library(dplyr)
library(patchwork)
##Example files can be downloaded via https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
#------------------------------User entries---------------------------------------------------------
#Set working directory for experiment
setwd("/Users/cburudpakdee/Desktop")
#Set path to features file
features_path <- "/Users/cburudpakdee/Desktop/GSM5931337_sc5_D18_BCMA_d0_features.tsv.gz"
#Set path to barcode file
barcode_path <- "/Users/cburudpakdee/Desktop/GSM5931337_sc5_D18_BCMA_d0_barcodes.tsv.gz"
#Set path to matrix file
matrix_path <- "/Users/cburudpakdee/Desktop/GSM5931337_sc5_D18_BCMA_d0_matrix.mtx.gz"
#Set minimum number of cells per gene
min_cells <- 50
#Set minimum number of features per cell
min_features <- 200
#Set number of top genes to compare
top_genes <- 2000
#-------------------------------End User entries and Run--------------------------------
cat("Setting seed")
set.seed(123)
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
cat("Statistics for cluster 1")
cluster1 <- head(FindMarkers(CAR_T,ident.1 = 1,min.pct = 0.25),n=100)
cluster1 <- rownames(cluster1)
cluster1<- data.frame(cluster1)
cluster1
cat("Statistics for cluster 2")
cluster2 <- head(FindMarkers(CAR_T,ident.1 = 2,min.pct = 0.25),n=100)
cluster2 <- row.names(cluster2)
cluster2<- data.frame(cluster2)
cluster2
gc()
cat("Statistics for cluster 3")
cluster3 <- head(FindMarkers(CAR_T,ident.1 = 3,min.pct = 0.25),n=100)
cluster3 <- row.names(cluster3)
cluster3 <- data.frame(cluster3)
cluster3

cat("Statistics for cluster 4")
cluster4 <- head(FindMarkers(CAR_T,ident.1 = 4,min.pct = 0.25),n=100)
cluster4 <- row.names(cluster4)
cluster4 <- data.frame(cluster4)
cluster4
gc()
cat("Statistics for cluster 5")
cluster5 <- head(FindMarkers(CAR_T,ident.1 = 5,min.pct = 0.25),n=100)
cluster5 <- row.names(cluster5)
cluster5 <- data.frame(cluster5)
cluster5
gc()
cat("Statistics for cluster 6")
cluster6 <- head(FindMarkers(CAR_T,ident.1 = 6,min.pct = 0.25),n=100)
cluster6 <- row.names(cluster6)
cluster6 <- data.frame(cluster6)
cluster6
gc()
cat("Statistics for cluster 7")
cluster7 <- head(FindMarkers(CAR_T,ident.1 = 7,min.pct = 0.25),n=100)
cluster7 <- row.names(cluster7)
cluster7 <- data.frame(cluster7)
cluster7
gc()
cat("Statistics for cluster 8")
cluster8 <- head(FindMarkers(CAR_T,ident.1 = 8,min.pct = 0.25),n=100)
cluster8 <- row.names(cluster8)
cluster8 <- data.frame(cluster8)
cluster8
gc()
cluster9 <- head(FindMarkers(CAR_T,ident.1 = 9,min.pct = 0.25),n=100)
cluster9 <- row.names(cluster9)
cluster9 <- data.frame(cluster9)
cluster9
gc()
cat("Statistics for cluster 10")
cluster10 <- head(FindMarkers(CAR_T,ident.1 = 10,min.pct = 0.25),n=100)
cluster10 <- row.names(cluster10)
cluster10 <- data.frame(cluster10)
cluster10
gc()
cat("Statistics for cluster 11")
cluster11 <- head(FindMarkers(CAR_T,ident.1 = 11,min.pct = 0.25),n=100)
cluster11 <- row.names(cluster11)
cluster11 <- data.frame(cluster11)
cluster11
gc()
cat("Statistics for cluster 12")
cluster12 <- head(FindMarkers(CAR_T,ident.1 = 12,min.pct = 0.25),n=100)
cluster12 <- row.names(cluster12)
cluster12 <- data.frame(cluster12)
cluster12
gc()
cat("Statistics for cluster 13")
cluster13 <- head(FindMarkers(CAR_T,ident.1 = 13,min.pct = 0.25),n=100)
cluster13 <- row.names(cluster13)
cluster13 <- data.frame(cluster13)
cluster13
gc()
cat("Statistics for cluster 14")
cluster14 <- head(FindMarkers(CAR_T,ident.1 = 14,min.pct = 0.25),n=100)
cluster14 <- row.names(cluster14)
cluster14 <- data.frame(cluster14)
cluster14
gc()


#markers <- FindAllMarkers(CAR_T,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
#markers %>% 
#  group_by(cluster) %>%
#  slice_max(n=2,order_by=avg_log2FC)
break
cat("Select and enter genes of interest")
cat("Enter in genes of interest after scanning list")
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","GZMB","GATA3","CCR7","PRF1","TIGIT"))
VlnPlot(CAR_T,features = c("CD4","CD8A","LAG3","CCR7","TNF","HLA-DRB1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","CTLA4","HLA-DRA","TNFRSF9","TNFRSF4","CD74"))
#trying iterating all genes... blah
#PCA1
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","H3C7","H3C8","H2BC7","H3C2","H2AC11"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","H2AC14","H3C10","CCL1","CCL3","H1-5"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","H1-2","H3C11","H4C3","H1-4","CXCL10"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","IL22","H2AC12","IFIT1","GNLY","PMCH"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","MT-ND6","CYBA","MT-CYB","MT-ND5","MT-ND2"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","MTDH","RBBP6","RPS2","ANKRD12","MT-ND1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","BTG1","SFPQ","ITGB7","PTPN7","PABPN1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","FYB1","TAP1","MT-ND4","MT-ATP8","SPOCK2"))
#PCA2
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","MKI67","CENPF","UBE2C","RRM2","GTSE1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","TOP2A","CDCA3","KIFC1","CDK1","BIRC5"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","ASPM","MXD3","PKMYT1","STMN1","AURKB"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","H2AX","CENPE","CDKN3","CKS1B","CCNA2"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","CD27","ISG20","FKBP11","GZMM","ERN1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","KDM5B","MZB1","ITGB7","CD7","CD40LG"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","MYDGF","CORO1B","PBXIP1","PRF1","LAT"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","TYMP","SEC11C", "SLC2A3","BTG1","HSPA5"))
#PCA3
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","CDC20","PLK1","TROAP","CENPE","CCNB1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","ASPM","KIF20A","CCNB2","KIF14","CENPF"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","PRR11","CENPA","PTTG1","PSRC1","ARL6IP1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","DLGAP5","KNSTRN","HMMR","PIMREG","BIRC5"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","GINS2","MCM7","MCM5","CLSPN","MCM3"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","FEN1","UHRF1","MCM2","PCNA","CDC45"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","CDCA7","MCM6","UNG","CDC6","TK1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","HELLS","CDK4","MCM4","DHFR","DUT"))
#PCA4
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","TCF7","PDE7B","SOX4","CCR7","GNG8"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","SPINT2","TNFRSF4","SESN3","NR3C1","BIRC3"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","C20orf204","PTGIR","HLA-DQA2","HLA-DQA1","RGS1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","DHRS3","TNFAIP3","STAG3","COL6A2","NFKB2"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","HOPX","NKG7","AQP3","ITGB7","MOSPD3"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","CDC25B","GZMA","ANXA2","CXCR3","SAMD3"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","TSPAN32","RAMP1","CTSC","PLEC","LGALS3"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","ZFP36","PLP2","NCR3","CST7","LGALS1"))
#PCA5
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","CCL5","CXCR6","PHLDA1","GZMB","IL13"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","ANXA1","PHLDA2","CSF2","CST7","TNFSF14"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","TNFRSF18","ATP8B4","NDFIP2","GADD45G","CCL3"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","XCL1","MAP3K8","CXCR3","CKLF","CCL4"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","FAM13A","AQP3","CCR7","TRABD2A","TCF7"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","CD27","NOSIP","TCEA3","S1PR1","SELL"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","GIMAP7","PLAC8","FHIT","KLF2","RAMP1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","ARMH1","MYC","GBP1","C1orf162","HAPLN3"))
#PCA6
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","IL9R","LIF","IRX3","NR4A1","CDKN2A"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","PHLDA2","SNCA","MMP25","CD33","TMEM273"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","HPGDS","ARL4A","MZB1","MBOAT7","PLPP1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","PDE7B","RGS16","OSM","GATA3","PTGDR2"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","CD8A","CD8B","CXCR3","HLA-DRB1","HLA-DRA"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","XCL1","HLA-DQA1","HLA-DRB5","HLA-DPA1","HLA-DQA2"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","NKG7","CCL5","LAG3","HOPX","HLA-DRB1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","HLA-DMA","HLA-DQB1","GZMB","TCF7","NELL2"))
#PCA7
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","AK4","BNIP3","CD8B","CD8A","BNIP3L"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","FAM162A","ALDOC","P4HA1","SLC16A3","GZMA"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","ENO2","PFKFB4","CTSW","NKG7","CD27"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","GPI","PLAC8","LAYN","ERO1A","XCL1"))
VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","HLA-DQA2","HLA-DQA1","HLA-DRA","HLA-DRB5","HLA-DPA1"))





VlnPlot(CAR_T,features = c("CAR-pCCL-BCMA","","","","",""))
print(CAr_T[["pca"]], dim=1:10,nFeatures=10)
cat("PC_ 1 
Positive:  H3C7, H3C8, H2BC7, H3C2, H2AC11, H2AC14, H3C10, CCL1, CCL3, H1-5 
	   H1-2, H3C11, H4C3, H1-4, CXCL10, IL22, H2AC12, IFIT1, GNLY, PMCH 
Negative:  MT-ND6, CYBA, MT-CYB, MT-ND5, MT-ND2, MTDH, RBBP6, RPS2, ANKRD12, MT-ND1 
	   BTG1, SFPQ, ITGB7, PTPN7, PABPN1, FYB1, TAP1, MT-ND4, MT-ATP8, SPOCK2 
PC_ 2 
Positive:  MKI67, CENPF, UBE2C, RRM2, GTSE1, TOP2A, CDCA3, KIFC1, CDK1, BIRC5 
	   ASPM, MXD3, PKMYT1, STMN1, AURKB, H2AX, CENPE, CDKN3, CKS1B, CCNA2 
Negative:  CD27, ISG20, FKBP11, GZMM, ERN1, KDM5B, MZB1, ITGB7, CD7, CD40LG 
	   MYDGF, CORO1B, PBXIP1, PRF1, LAT, TYMP, SEC11C, SLC2A3, BTG1, HSPA5 
PC_ 3 
Positive:  CDC20, PLK1, TROAP, CENPE, CCNB1, ASPM, KIF20A, CCNB2, KIF14, CENPF 
	   PRR11, CENPA, PTTG1, PSRC1, ARL6IP1, DLGAP5, KNSTRN, HMMR, PIMREG, BIRC5 
Negative:  GINS2, MCM7, MCM5, CLSPN, MCM3, FEN1, UHRF1, MCM2, PCNA, CDC45 
	   CDCA7, MCM6, UNG, CDC6, TK1, HELLS, CDK4, MCM4, DHFR, DUT 
PC_ 4 
Positive:  TCF7, PDE7B, SOX4, CCR7, GNG8, SPINT2, TNFRSF4, SESN3, NR3C1, BIRC3 
	   C20orf204, PTGIR, HLA-DQA2, HLA-DQA1, RGS1, DHRS3, TNFAIP3, STAG3, COL6A2, NFKB2 
Negative:  HOPX, NKG7, AQP3, ITGB7, MOSPD3, CDC25B, GZMA, ANXA2, CXCR3, SAMD3 
	   TSPAN32, RAMP1, CTSC, PLEC, LGALS3, ZFP36, PLP2, NCR3, CST7, LGALS1 
PC_ 5 
Positive:  CCL5, CXCR6, PHLDA1, GZMB, IL13, ANXA1, PHLDA2, CSF2, CST7, TNFSF14 
	   TNFRSF18, ATP8B4, NDFIP2, GADD45G, CCL3, XCL1, MAP3K8, CXCR3, CKLF, CCL4 
Negative:  FAM13A, AQP3, CCR7, TRABD2A, TCF7, CD27, NOSIP, TCEA3, S1PR1, SELL 
	   GIMAP7, PLAC8, FHIT, KLF2, RAMP1, ARMH1, MYC, GBP1, C1orf162, HAPLN3 
PC_ 6 
Positive:  IL9R, LIF, IRX3, NR4A1, CDKN2A, PHLDA2, SNCA, MMP25, CD33, TMEM273 
	   HPGDS, ARL4A, MZB1, MBOAT7, PLPP1, PDE7B, RGS16, OSM, GATA3, PTGDR2 
Negative:  CD8A, CD8B, CXCR3, HLA-DPB1, HLA-DRA, XCL1, HLA-DQA1, HLA-DRB5, HLA-DPA1, HLA-DQA2 
	   NKG7, CCL5, LAG3, HOPX, HLA-DRB1, HLA-DMA, HLA-DQB1, GZMB, TCF7, NELL2 
PC_ 7 
Positive:  AK4, BNIP3, CD8B, CD8A, BNIP3L, FAM162A, ALDOC, P4HA1, SLC16A3, GZMA 
	   ENO2, PFKFB4, CTSW, NKG7, CD27, GPI, PLAC8, LAYN, ERO1A, XCL1 
Negative:  HLA-DQA2, HLA-DQA1, HLA-DRA, HLA-DRB5, HLA-DPA1, HLA-DPB1, HLA-DQB1, SELL, PRDX1, CTSH 
	   SEC11C, HLA-DRB1, NOP16, SDF2L1, AQP3, CD28, DDX21, S1PR1, MRTO4, IFRD2 
PC_ 8 
Positive:  CD8B, CD8A, CTSW, TRIM22, LRRN3, XCL1, FCER1G, FDXR, TCF7, PHPT1 
	   PRF1, NKG7, CD7, AGTRAP, IL9R, BAX, RPS27L, MT2A, TNFSF8, CDKN2A 
Negative:  AK4, BNIP3, FAM162A, P4HA1, ALDOC, SLC16A3, BNIP3L, GPI, PFKFB4, ENO2 
	   HLA-DRB1, ERO1A, HLA-DRB5, GAPDH, NOL3, ENO1, LDHA, HBEGF, PGK1, IL2RA 
PC_ 9 
Positive:  MTRNR2L12, MTRNR2L8, PBXIP1, MTRNR2L1, TNFSF8, CD28, H2AC21, H1-4, NUMA1, RNF213 
	   H2AC4, CDKN1A, MYH9, ESCO2, PKMYT1, FDXR, HBEGF, H1-5, SPTAN1, AHNAK 
Negative:  CD8A, CD8B, FABP5, CDC20, NOP16, XCL1, CTSW, MRTO4, PAICS, AK4 
	   FAM162A, CCNB1, LDHA, DCTPP1, MYC, MIF, FAM216A, DDX21, HSPE1, EBNA1BP2 
PC_ 10 
Positive:  ANXA1, FDXR, RPS27L, DDB2, TNFRSF18, HSPB1, NDFIP2, ACP5, SELENOW, ITM2A 
	   ALOX5AP, TNFSF8, PIM2, BAX, KLRB1, CTSC, TIMP1, RAD51C, TSC22D3, SRGN 
Negative:  CD8A, CD8B, MYH9, ADAM19, MTRNR2L8, MTRNR2L12, NR4A1, SGK1, LIF, PRKX 
	   FLNA, H1-4, XCL1, MTRNR2L1, MACF1, HBEGF, PPP1R16B, NFKB2, SF3B3, HMGCS1 
")

#
FeaturePlot(CAR_T,features=c("CCR7","GZMB","IL7R","CAR-pCCL-BCMA","CD3D","CD27"),reduction= "tsne")
FeaturePlot(CAR_T,features=c("CD4","CD8A","LAG3","CCR7","TNF","HLA-DRB1"),reduction= "tsne")
FeaturePlot(CAR_T,features=c("CCR7","GZMB","IL7R","CAR-pCCL-BCMA","CD3D","CD27"),reduction= "umap")
FeaturePlot(CAR_T,features=c("CD4","CD8A","LAG3","CCR7","TNF","HLA-DRB1"),reduction= "umap")            
FeaturePlot(CAR_T,features=c("CCR7","GZMB","IL7R","CAR-pCCL-BCMA","CD3D","CD27"),reduction= "pca")
FeaturePlot(CAR_T,features=c("CD4","CD8A","LAG3","CCR7","TNF","HLA-DRB1"),reduction= "pca")            
cat("Analysis is complete")
cluster_ids <- c("noise","DNA_replication","cytolytic_granule/allograft_rejection","T_cell_receptor_pathway","ribosomal translation","cell cycle/centromere","innate immune response",
                 "isopeptide bonds,blot clots","____________actin filament binding","microtubule binding","TNF/NF-kappaB pathway","Epstein Barr infection, MHC class1, MHC class 2","Innate response. Response to virus","mitochondrion","mitochondrion")
names(cluster_ids) <- levels(CAR_T)
CART <- RenameIdents(CAR_T,cluster_ids)
DimPlot(CART,reduction = "tsne",label = TRUE, pt.size = 0.05)+ NoLegend()
#---------------------------------------------------------------------------
