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
counts <- readMM("matrix.mtx.gz")
genes <- read_tsv("genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
head(gene_ids)
cell_ids <- read_tsv("barcodes.tsv.gz", col_names = FALSE)$X1
head(cell_ids)
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
matrix<-data.frame(counts)
names<-c(paste0("norm_", 1:11883))
colnames(matrix)<-names
#-------------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/gbm/R2_T')
counts2 <- readMM("matrix.mtx.gz")
genes2 <- read_tsv("genes.tsv.gz", col_names = FALSE)
gene_ids2 <- genes2$X2
head(gene_ids2)
cell_ids2 <- read_tsv("barcodes.tsv.gz", col_names = FALSE)$X1
head(cell_ids2)
rownames(counts2) <- gene_ids2
colnames(counts2) <- cell_ids2
matrix2<-data.frame(counts2)
names2<-c(paste0("tumor_", 1:15094))
colnames(matrix2)<-names2
#-----Even out files here
matrix<-matrix[,c(1:1000)]
matrix2<-matrix2[,c(1:1000)]
#-------------------------------------------------------------------------------
Matrix<-cbind(matrix,matrix2)
#---must add 1 to data frame
Matrix<-Matrix+1
rm(counts,counts2,genes,genes2,matrix,matrix2)
gc()
names<-colnames(Matrix)
condition<-c(replicate(1000,'normal'),replicate(1000,'tumor'))
str(condition)
dim(Matrix)
type<-replicate(2000, 'paired')
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
deseq<-estimateSizeFactors(deseq)
DE<-DESeq(deseq,
          test = c('LRT'),
          useT = TRUE,
          fitType = 'glmGamPoi',
          minmu = 1e-6,
          minReplicatesForReplace = Inf,
          reduced = ~ 1)
resLRT <- results(DE)
resLRT<-data.frame(resLRT)
resLRT
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results_df<-data.frame(resLRT)
results_df<-results_df[order(results_df$padj, decreasing=FALSE),]
results_df<-subset(results_df,padj< 0.05)
summary(results_df)
#-------------------------------------------------------------------------------
EnhancedVolcano(resLRT,
                lab = row.names(resLRT),
                x='log2FoldChange',
                y='pvalue',
                title = 'Differential gene expressions of glioblastoma tumors compared with adjacent tissues',
                subtitle = 'Likelihood ratio test | n=2000 cells',
                legendLabels = NULL,
                legendIconSize = -1,
                legendPosition = 'bottom',
                pCutoff = 0.001,
                FCcutoff = 0.2,
                shape = 1,
                ylim = c(0,80),
                xlim = c(-1,1))

