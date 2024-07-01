library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
library(Matrix)
library(DESeq2)
library(EnhancedVolcano)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/gbm/R2_N')
# Read in `matrix.mtx`
counts <- readMM("matrix.mtx.gz")
# Read in `genes.tsv`
genes <- read_tsv("genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
head(gene_ids)
# Read in `barcodes.tsv`
cell_ids <- read_tsv("barcodes.tsv.gz", col_names = FALSE)$X1
head(cell_ids)
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
matrix<-data.frame(counts)
names<-c(paste0("norm_", 1:11883))
colnames(matrix)<-names
#-------------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/gbm/R2_T')
# Read in `matrix.mtx`
counts2 <- readMM("matrix.mtx.gz")
# Read in `genes.tsv`
genes2 <- read_tsv("genes.tsv.gz", col_names = FALSE)
gene_ids2 <- genes2$X2
head(gene_ids2)
# Read in `barcodes.tsv`
cell_ids2 <- read_tsv("barcodes.tsv.gz", col_names = FALSE)$X1
head(cell_ids2)
rownames(counts2) <- gene_ids2
colnames(counts2) <- cell_ids2
matrix2<-data.frame(counts2)
names2<-c(paste0("tumor_", 1:15094))
colnames(matrix2)<-names2
#-----Even out files here
#matrix2<-matrix2[,c(1:11883)]
matrix<-matrix[,c(1:2000)]
matrix2<-matrix2[,c(1:2000)]
#-------------------------------------------------------------------------------
Matrix<-cbind(matrix,matrix2)
Matrix<-Matrix+1
#rm(counts,counts2,genes,genes2,matrix,matrix2)
#gc()
break
names<-colnames(Matrix)
dim(matrix)
dim(matrix2)
condition<-c(replicate(2000,'normal'),replicate(2000,'tumor'))
str(condition)
dim(Matrix)
type<-replicate(4000, 'paired')
#rm(counts,counts2,genes,genes2,matrix,matrix2)
gc()
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
head(coldata)
colnames(Matrix)[1:5]
row.names(coldata)[1:5]
gc()
#-------------------------------------------------------------------------------
deseq<-DESeqDataSetFromMatrix(countData = Matrix,
                              colData = coldata,
                              design = ~ condition)
gc()
deseq
#rm(Matrix)
#gc()
DE<-DESeq(deseq,
          test = c('Wald'),
          parallel = TRUE)
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
results_df<-data.frame(results)
results_df<-results_df[order(results_df$padj, decreasing=FALSE),]
results_df<-subset(results_df,padj< 0.05)
summary(results_df)
#-------------------------------------------------------------------------------
EnhancedVolcano(results,
                lab = row.names(results),
                x='log2FoldChange',
                y='pvalue',
                title = "",
                subtitle = '',
                legendLabels = NULL,
                legendIconSize = -1,
                legendPosition = 'bottom',
                pCutoff = 0.001,
                FCcutoff = 0.2,
                shape = 1,
                ylim = c(0,50),
                xlim = c(-0.5,0.5))

