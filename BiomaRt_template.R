library(biomaRt)
library(dplyr)
setwd('/home/em_b/Desktop/FCCC/t3_edits')
matrix <- read.csv('quant_hg38.csv')
matrix<-matrix[,-c(2:3)]
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
dim(matrix)
geneid <- matrix$Name
head(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_transcript_id_version','external_gene_name','entrezgene_id'),
              filters = 'ensembl_transcript_id_version',
              values = geneid,
              mart = ensembl)
head(genes)
df <- matrix[matrix$Name %in% unique(genes$ensembl_transcript_id_version),]
dim(df)
dim(matrix)
genes$ensembl_transcript_id_version
colnames(genes)[1]<- 'Gene'
colnames(df)[1]<-'Gene'
Ftable <- merge(genes, df, by="Gene")
Ftable<-Ftable[order(Ftable$external_gene_name),]
colnames(Ftable)[1] <- c('transcript_version')
colnames(Ftable)[2] <- c('Gene')
rownames(Ftable)<-make.names(Ftable$Gene,unique = TRUE)
Ftable<-cbind(rownames(Ftable),Ftable)
colnames(Ftable)[1] <- c('tracking_id')
head(Ftable)
#matrix<-data.frame(Ftable)
#------Remove duplicate genes by p_value
#matrix<-matrix%>%
#  group_by(Gene)%>%
#  arrange(p_value)%>%
#  slice(1)
#matrix<-data.frame(matrix)
#str(matrix)
break
write.csv(Ftable,file = 'biomaRt.csv')
#----if using older mm10 mouse assembly
ensembl<-useEnsembl(biomart = 'genes',
                    dataset = 'mmusculus_gene_ensembl',
                    version=102)
