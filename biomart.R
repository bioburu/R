library(biomaRt)
library(dplyr)
setwd('/home/AST_OPC_files')
matrix <- read.csv('salmon_out.csv')
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
dim(matrix)
geneid <- matrix$Gene
head(geneid)
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
colnames(df)[1]<-'Gene'
Ftable <- merge(genes, df, by="Gene")
Ftable<-Ftable[order(Ftable$external_gene_name),]
ensembl<-useEnsembl(biomart = 'genes',
                    dataset = 'mmusculus_gene_ensembl',
                    version=102)
