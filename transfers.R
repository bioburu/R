#-----get RNA
features_path <- 'GSM5952326_MB3076_genes.tsv.gz'
barcodes_path <- 'GSM5952326_MB3076_barcodes.tsv.gz'
matrix_path <- 'GSM5952326_MB3076_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
rna <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'rna')
table(rna@active.ident)
table(atac@active.ident)
atac<-subset(x = atac, downsample = 1370)
#--------------------
gc()
break 
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 11000 & percent.mt <15)
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
rna <- NormalizeData(rna)
gc()
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
VariableFeatures(rna)
top1000 <- head(VariableFeatures(rna), 1000)
top1000
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
gc()
p1<-DimPlot(rna, reduction = "umap",dims = c(1,2),cols = 'red')+ggtitle('scRNAseq')
p2<-DimPlot(atac, reduction = "umap.atac",dims = c(1,2),cols = 'black')+ggtitle('scATACseq')
p1+p2
break
