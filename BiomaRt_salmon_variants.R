library(biomaRt)
library(dplyr)
setwd('/home/em_b/Downloads')
matrix <- read.csv('matrix.csv')
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
dim(matrix)
geneid <- matrix$Gene
head(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_transcript_id_version','external_gene_name','entrezgene_id'),
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
colnames(Ftable)[1] <- c('transcript_version')
colnames(Ftable)[2] <- c('Gene')
rownames(Ftable)<-make.names(Ftable$Gene,unique = TRUE)
Ftable<-cbind(rownames(Ftable),Ftable)
colnames(Ftable)[1] <- c('tracking_id')
head(Ftable)
matrix<-Ftable
library(dplyr)
str(matrix)
matrix<-matrix %>%
  group_by(Gene) %>%
  summarize(T3_1 = sum(T3_1),
            T3_2 = sum(T3_2),
            T3_3 = sum(T3_3),
            PBS_1 = sum(PBS_1),
            PBS_2 = sum(PBS_2),
            PBS_3 = sum(PBS_3))
break 
write.csv(matrix,file = 'matrix.csv')
#----if using older mm10 mouse assembly
ensembl<-useEnsembl(biomart = 'genes',
                    dataset = 'mmusculus_gene_ensembl',
                    version=102)
