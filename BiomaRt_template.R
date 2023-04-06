library(biomaRt)
setwd("/home/amp_prog/Desktop")
matrix <- read.delim('screening_1.txt')
matrix<-matrix[,-c(2:4)]
colnames(matrix)<-c('Gene','relapse_2')
#---------------------------------------------------------------------
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
genes <-getBM(attributes = c('ensembl_transcript_id_version','external_gene_name'),
      filters = 'ensembl_transcript_id_version',
      values = geneid,
      mart = ensembl)
head(genes)
df <- matrix[matrix$Gene %in% unique(genes$ensembl_transcript_id_version),]
dim(df)
dim(matrix)
genes$ensembl_transcript_id_version
colnames(genes)[1]<- 'Gene'
Ftable <- merge(genes, df, by="Gene")
Ftable<- Ftable[,-c(1)]
colnames(Ftable)[1] <- c('Gene')
View(Ftable)
write.csv(Ftable,file = 'relapse_2.csv')

