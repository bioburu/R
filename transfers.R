
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
setwd('/Users/burudpc/Desktop/GSE198565_RAW')
matrix<-Matrix::readMM('GSM5952324_MB3076_ATAC_matrix.mtx.gz')
barcodes <-readLines('GSM5952324_MB3076_ATAC_barcodes.tsv.gz')
peaks<-read.table('GSM5952324_MB3076_ATAC_peaks.bed.gz',sep='\t')
#--concatenate all columns into list of strings for rownames
peaknames<-paste(peaks$V1,peaks$V2,peaks$V3,sep = '-')
colnames(matrix)<-barcodes
row.names(matrix)<-peaknames
head(matrix)
atac <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'atac')
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
head(atac)
DimPlot(atac,label = TRUE)+ggtitle('scATACseq')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
head(annotations)
#--change chromosome annotation style 
seqlevelsStyle(annotations) <- "UCSC"
#---reassign genome annotations if necessary
genome(annotations) <- "hg38"
genome(annotations)

atac
df<-data.frame(atac[['RNA']]$counts)

break 
pbmc.atac@meta.data
#---add annotations to atac file
Annotation(pbmc.atac) <- annotations
#chrom_assay<-CreateChromatinAssay(counts = matrix,sep = c(':','-'))
break
summary(x@active.ident)
x<-subset(x = x, downsample = 3160)
summary(x@active.ident)
#---convert into Assay5 obj
#pbmc.rna
#pbmc.rna[['RNA']]<-as(pbmc.rna[['RNA']],Class = 'Assay5')
#head(pbmc.rna@meta.data)
#---filter out 'filtered' category
#pbmc.rna<-subset(pbmc.rna,seurat_annotations!='filtered')
#head(pbmc.rna@meta.data)
#pbmc.rna <- NormalizeData(pbmc.rna)
#pbmc.rna <- FindVariableFeatures(pbmc.rna)
#pbmc.rna <- ScaleData(pbmc.rna)
#pbmc.rna <- RunPCA(pbmc.rna)
#pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)
gc()
#VlnPlot(pbmc.rna,features = c('PTPRC','CD3D'))
#DimPlot(pbmc.rna,group.by = 'seurat_annotations',label = TRUE)+ggtitle('scRNAseq')
#-------------------------------------------------------------------------------
load(file = '/home/deviancedev/Desktop/scATACseq/pbmcMultiome.SeuratData/data/pbmc.atac.rda')
head(pbmc.atac@meta.data)
#---remove 'filtered' category from annotations
pbmc.atac<-subset(pbmc.atac,seurat_annotations!='filtered')
head(pbmc.atac@meta.data)
# We exclude the first dimension as this is typically correlated with sequencing depth
#---normalized by term-frequency inverse-document-frequency
pbmc.atac <- RunTFIDF(pbmc.atac)
#----find top features based on number counts/feature
#---q0 denotes top 100% of most common features
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
#---apply Singular Value Decomposion to matrix
pbmc.atac <- RunSVD(pbmc.atac)
#---umap 
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
head(pbmc.atac)
break 
#DimPlot(pbmc.atac,group.by = 'seurat_annotations',label = TRUE)+ggtitle('scATACseq')
DimPlot(pbmc.atac,label = TRUE)+ggtitle('scATACseq')
# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
head(annotations)
#--change chromosome annotation style 
seqlevelsStyle(annotations) <- "UCSC"
#---reassign genome annotations if necessary
genome(annotations) <- "hg38"
genome(annotations)
head(data.frame(pbmc.atac[['ATAC']]$counts))

break 
pbmc.atac@meta.data
#---add annotations to atac file
Annotation(pbmc.atac) <- annotations
p1 <- DimPlot(pbmc.rna, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2
# get most variable gene activity from the RNAseq file of the ATACseq file
gene.activities <- GeneActivity(pbmc.atac, features = VariableFeatures(pbmc.rna))
gc()
View(data.frame(gene.activities))
Annotation(pbmc.atac)
# add gene activities as a new assay
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# normalize gene activities
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac@assays$ACTIVITY
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
gc()
#----predict subsets of ATAC data using anchors from RNA data
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$seurat_annotations,
                                     weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)
View(celltype.predictions)
#---add cell type predictions to ATAC metadata
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
head(pbmc.atac@meta.data)

break 
#---add annotation_correct=TRUE/FALSE to ATAC metadata
#---only works if using Multiome kits
pbmc.atac$annotation_correct <- pbmc.atac$predicted.id == pbmc.atac$seurat_annotations
head(pbmc.atac@meta.data)


p1 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(pbmc.atac, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2

predictions <- table(pbmc.atac$seurat_annotations, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p1


break 
correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + 
  theme_cowplot() + 
  scale_fill_discrete(name = "Annotation Correct",labels = c(paste0("FALSE (n = ", incorrect, ")"),paste0("TRUE (n = ", correct, ")"))) + 
  scale_color_discrete(name = "Annotation Correct",labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + 
  xlab("Prediction Score")
p1 + p2
#
View(data.frame(pbmc.atac[['ATAC']]$counts))
break 
predictions <- table(pbmc.atac$seurat_annotations, pbmc.atac$predicted.id)
predictions
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions
predictions <- as.data.frame(predictions)
predictions
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                                                                                                                                  labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2
