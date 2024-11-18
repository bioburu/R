#remotes::install_github("mikelove/fishpond")
library(fishpond)
library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
library(Seurat)
library(Matrix)
library(readr)
files<-'/home/mm10_quant'
x<-loadFry(files,outputFormat = 'scRNA',nonzero = FALSE,quiet = FALSE)
x
gc()
counts <- assay(x, "counts")
counts
dim(counts)
df<-data.frame(counts)
gc()
matrix<-df
matrix<-cbind(row.names(matrix),matrix)
colnames(matrix)[1]<-'Gene'
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl<-useEnsembl(biomart = 'genes',
                   dataset = 'mmusculus_gene_ensembl',
                   version=102)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
dim(matrix)
geneid <- matrix$Gene
head(geneid)
dim(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_gene_id','external_gene_name','chromosome_name','p_value','start_position','end_position','description'),
              filters = 'ensembl_gene_id',
              values = geneid,
              mart = ensembl)
head(genes)
genes<-genes[!duplicated(genes$ensembl_gene_id), ]
df <- matrix[matrix$Gene %in% unique(genes$ensembl_gene_id),]
dim(df)
dim(matrix)
colnames(genes)[1]<- 'Gene'
Ftable <- merge(genes, df, by="Gene")
gc()
colnames(Ftable)[1] <- c('ensembl_id')
rm(attributes,counts,df,ensembl,filters,genes,matrix,x)
gc()
Ftable<-Ftable[!duplicated(Ftable$external_gene_name), ]
#write.csv(Ftable,file = 'GSM8004749.Sciatic.nerve.CD45p.spon.autoimmune.neuropathy.csv')
features<-Ftable[,c(1:2)]
matrix.mtx<-Ftable[,-c(1:7)]
matrix.mtx<-as(matrix.mtx,'sparseMatrix')
barcodes<-data.frame(colnames(Ftable)[8:8800])
setwd('/home/drive')
write_tsv(features,file='features.tsv',col_names = NA)
write_tsv(barcodes,file='barcodes.tsv',col_names = NA)
writeMM(matrix.mtx,file='matrix.mtx')
