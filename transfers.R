library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
setwd('/Users/burudpc/Desktop/GSE198565_RAW')
#----get ATAC
matrix<-Matrix::readMM('GSM5952324_MB3076_ATAC_matrix.mtx.gz')
barcodes <-readLines('GSM5952324_MB3076_ATAC_barcodes.tsv.gz')
peaks<-read.table('GSM5952324_MB3076_ATAC_peaks.bed.gz',sep='\t')
#--concatenate all columns into list of strings for rownames
peaknames<-paste(peaks$V1,peaks$V2,peaks$V3,sep = '-')
colnames(matrix)<-barcodes
row.names(matrix)<-peaknames
head(matrix)
#---------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
head(annotations)
#--change chromosome annotation style 
seqlevelsStyle(annotations) <- "UCSC"
#---reassign genome annotations if necessary
genome(annotations) <- "hg38"
genome(annotations)
head(annotations)
atac<-CreateChromatinAssay(counts = matrix,sep = c(':','-'),annotation = annotations)
atac <- CreateSeuratObject(counts=atac,min.cells=20,min.features=200,project = 'atac')
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
head(atac)
DimPlot(atac,label = TRUE)+ggtitle('scATACseq')
dim(atac)
