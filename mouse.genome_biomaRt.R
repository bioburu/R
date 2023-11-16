library(biomaRt)
library(devtools)
#devtools::install_version("dbplyr", version = "2.3.4")
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/alignments/finals')
matrix <- read.csv('master_ensembl.csv')
#----If gene or transcript versions are not needed use:
matrix$Gene
matrix$Gene<-sub('\\..*','',matrix$Gene)
matrix$Gene
geneid <- matrix$Gene
head(geneid)
#---------------------------------------------------------------------
head(listMarts(host='https://www.ensembl.org'),10)
head(biomaRt::listDatasets(biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")), 10)     
head(biomaRt::listAttributes(biomaRt::useDataset(dataset = "mmusculus_gene_ensembl", mart= useMart("ENSEMBL_MART_ENSEMBL",host= "https://www.ensembl.org"))), 10)
mart <- biomaRt::useDataset(dataset = "mmusculus_gene_ensembl",mart= useMart("ENSEMBL_MART_ENSEMBL",host= "https://www.ensembl.org"))
genes <-getBM(attributes = c('ensembl_transcript_id','external_gene_name'),
              filters = 'ensembl_transcript_id',
              values = geneid,
              mart = mart)
#-----------------------------------------------------------------
df <- matrix[matrix$Gene %in% unique(genes$ensembl_transcript_id),]
dim(df)
dim(matrix)
genes$ensembl_transcript_id
colnames(genes)[1]<- 'Gene'
Ftable <- merge(genes, df, by="Gene")
Ftable<- Ftable[,-c(1)]
colnames(Ftable)[1] <- c('Gene')
break 
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/alignments/finals')
write.csv(Ftable,file = 'master_gene.list.csv')
