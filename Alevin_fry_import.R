#remotes::install_github("mikelove/fishpond")
library(fishpond)
library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
library(Seurat)
files<-'/home/deviancedev/Desktop/drive_jan2024/Alevin_Salmon/af_tutorial/alevin_pilot_quant'
x<-loadFry(files,outputFormat = 'scRNA',nonzero = FALSE,quiet = FALSE)
x
#----
gc()
counts <- assay(x, "counts")
counts
df<-data.frame(counts)
gc()
#------
matrix<-df
matrix<-cbind(row.names(matrix),matrix)
colnames(matrix)[1]<-'Gene'
#-----
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart = ensembl)
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
df <- matrix[matrix$Gene %in% unique(genes$ensembl_gene_id),]
dim(df)
dim(matrix)
colnames(genes)[1]<- 'Gene'
Ftable <- merge(genes, df, by="Gene")
gc()
colnames(Ftable)[1] <- c('ensembl_id')
break
setwd('/home/deviancedev/Desktop/drive_jan2024/Alevin_Salmon/af_tutorial')
write.csv(Ftable,file = 'alevin_pbmc_scRNAseq.csv')
