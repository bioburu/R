library(biomaRt)
library(dplyr)
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/RNAseq_exosomes/run_2')
matrix <- read.csv('exosomes_run2.csv')
matrix<-cbind(matrix$Gene,matrix)
matrix$`matrix$Gene`
matrix$`matrix$Gene`<-sub('\\..*','',matrix$`matrix$Gene`)
matrix$`matrix$Gene`
#---set host to https://useast.ensembl.org or https://uswest.ensembl.org
#---------------------------------------------------------------------
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
dim(matrix)
geneid <- matrix$`matrix$Gene`
head(geneid)
dim(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_transcript_id','external_gene_name','chromosome_name','p_value','start_position','end_position','description'),
              filters = 'ensembl_transcript_id',
              values = geneid,
              mart = ensembl)
head(genes)
df <- matrix[matrix$`matrix$Gene` %in% unique(genes$ensembl_transcript_id),]
dim(df)
dim(matrix)
genes$ensembl_transcript_id
colnames(genes)[1]<- 'matrix$Gene'
Ftable <- merge(genes, df, by="matrix$Gene")
Ftable<- Ftable[,-c(1)]
colnames(Ftable)[7] <- c('transcript_version')
colnames(Ftable)[1] <- c('Gene')
head(Ftable)
matrix<-data.frame(Ftable)
#------Remove duplicate genes by p_value
matrix<-matrix%>%
  group_by(Gene)%>%
  arrange(p_value)%>%
  slice(1)
matrix<-data.frame(matrix)
str(matrix)
break
write.csv(matrix,file = 'exosomes_rnaseq_run2.csv')
#----if using older mm10 mouse assembly
ensembl<-useEnsembl(biomart = 'genes',
                    dataset = 'mmusculus_gene_ensembl',
                    version=102)
