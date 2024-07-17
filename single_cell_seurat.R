library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
cat('GSE234832. Patient-derived brain metastasis models')
setwd('/home/em_b/work_stuff/brain_metastasis/GSE234832_RAW')
barcodes_path <- 'GSM7475327_LUBMET7.barcodes.tsv.gz'
features_path <- 'GSM7475327_LUBMET7.features.tsv.gz'
matrix_path <- 'GSM7475327_LUBMET7.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lung_cancer')
#---------------------------------------------------------------------------
setwd('/home/em_b/work_stuff/brain_metastasis/GSE234832_RAW')
barcodes_path <- 'GSM7475328_LUBMET1.barcodes.tsv.gz'
features_path <- 'GSM7475328_LUBMET1.features.tsv.gz'
matrix_path <- 'GSM7475328_LUBMET1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lung_cancer')
data<-merge(data,data2)
table(data@active.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <10)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
top1000
gc()
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
gc()
data <- ScaleData(data, features = all.genes)
gc()
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.6)
data <- RunUMAP(data, dims = 1:15)
data <-JoinLayers(data)
DimPlot(data,
        reduction = 'umap',
        label=TRUE,
        label.size = 6,
        repel = TRUE,
        label.box = FALSE,
        raster=FALSE,
        pt.size = 0.5,
        seed=1,
        cols.highlight = c('grey'),
        dims = c(1,2))
break 
FindMarkers(data,
            ident.1 = '7',
            ident.2 = '5',
            logfc.threshold=0,
            only.pos = TRUE,
            test.use='bimod')
#---tester
VlnPlot(data, features = c('ULK1', 'ATG13', 'FIP200', 'ATG101', 'ATG9A'),cols = c())

#-------------------------------------------------------------------------------
#---cluster 5 T-cells
VlnPlot(data, features = c('PTPRC','TRAC','TRBC1','CD3D','CD3E','CD8A','RORA'),cols = c())

#----cluster 2/8 Type II macro/microglial
VlnPlot(data,
        features = c('PTPRC','CD14','MRC1','MILR1','CD163',
                     'LYZ','CD86','CD68','FCGR3A','ITGAM',
                     'ITGAX','CSF1R'))

#----clusters 4,6,9 oligodendrocytes and astrocytes
VlnPlot(data,
        features = c('NCAM1','OLIG1','OLIG2','SOX8','SOX9','SOX10','FGFR2',
                     'FGF1','LAMP2','GFAP','S100B','NDRG2','SLC1A3'))

#---cluster 10 Fibroblasts
VlnPlot(data, features = c('COL1A1','COL6A1','COL1A2','COL4A1','FN1','IGFBP7','TIMP3'),cols = c())

#---cluster 0,1,3,7 
VlnPlot(data, features = c('EPCAM','CD24','SOX17'),cols = c())

#----cytokines
VlnPlot(data, features = c('IL1B','IL2RG','IL6R','IL6ST','CXCL8',
                           'IL10RA','IL18','IL32','TNF','TNFAIP3',
                           'TNFSF10','TGFB1','CSF1','CSF1R'),cols = c())
new_idents<-c('tumor','tumor','M2mac','tumor','OPC','T_cells','OPC','tumor','M1macro','OPC','FB')
names(new_idents) <- levels(data)
data <- RenameIdents(data, new_idents)
levels(data)
break 
#-------------------------------------------
markers<-FindMarkers(data,
            ident.1 = '7',
            ident.2 = '5',
            logfc.threshold=0,
            only.pos = TRUE,
            test.use='bimod')
#-----Define the clusters 
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
geneid <- row.names(markers)
head(geneid)
genes <-getBM(attributes = c('external_gene_name','entrezgene_id'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
go1<-data.frame(go)
go1<-go1[order(go1$Count, decreasing=TRUE),]
#-----------------------------------------------------------------------
GO_search<-select(org.Hs.eg.db,
                               keytype = 'GOALL',
                               keys = 'GO:0198738',
                               columns = c('SYMBOL','GENENAME','ENTREZID'))
list<-GO_search$SYMBOL
list
DoHeatmap(
  data,
  features = list,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)

#-----If subsetting use
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD45.groups')
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
#-------------CD45+ CD19+  ---------------------------------------------------------
data$CD19.groups <- 'CD19.pos'
data$CD19.groups[WhichCells(data, expression= CD19 < 0.1)] <- 'CD19.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD19.groups')
head(data@meta.data)
data <- subset(data, subset = CD19.groups != "CD19.neg")
gc()
#-------------CD45+ CD19+ TRAC- -------------------------------------------
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.1)] <- 'TRAC.neg'
DimPlot(data, reduction = 'pca',split.by = 'TRAC.groups')
head(data@meta.data)
data <- subset(data, subset = TRAC.groups != "TRAC.pos")
gc()
#-------------CD45+ CD19+ TRAC- CD14-     -------------------------------------------
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.1)] <- 'CD14.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD14.groups')
head(data@meta.data)
data <- subset(data, subset = CD14.groups != "CD14.pos")
gc()
#-------------CD45+ CD19+ TRAC- CD14- MKI67+    -------------------------------------------
data$MKI67.groups <- 'MKI67.pos'
data$MKI67.groups[WhichCells(data, expression= MKI67 < 0.1)] <- 'MKI67.neg'
DimPlot(data, reduction = 'pca',split.by = 'MKI67.groups')
head(data@meta.data)
data <- subset(data, subset = MKI67.groups != "MKI67.neg")
gc()
VlnPlot(data, features = c('PTPRC','CD19','TRAC','CD14','MKI67'),pt.size=0.1)

