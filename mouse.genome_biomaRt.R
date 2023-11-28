library(biomaRt)
library(devtools)
#---if doesnt work, reinstall dbplyr and restart Rstudio
#devtools::install_version("dbplyr", version = "2.3.4")
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/alignments/finals')
matrix <- read.csv('quant_out_unordered.csv')
#----If gene or transcript versions are not needed use:
matrix$Gene
matrix$Gene<-sub('\\..*','',matrix$Gene)
matrix$Gene
geneid <- matrix$Gene
#---set host to https://useast.ensembl.org or https://uswest.ensembl.org
#---------------------------------------------------------------------
#head(listMarts(host='https://www.ensembl.org'),10)
#head(biomaRt::listDatasets(biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")), 10)     
#head(biomaRt::listAttributes(biomaRt::useDataset(dataset = "mmusculus_gene_ensembl", mart= useMart("ENSEMBL_MART_ENSEMBL",host= "https://www.ensembl.org"))), 10)
mart <- biomaRt::useDataset(dataset = "mmusculus_gene_ensembl",mart= useMart("ENSEMBL_MART_ENSEMBL",host= "https://www.ensembl.org"))
listFilters(mart)
listAttributes(mart)
listDatasets(mart)
genes <-getBM(attributes = c('ensembl_transcript_id','external_gene_name','chromosome_name','start_position','end_position'),
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
Ftable<-Ftable[order(Ftable$Gene), ]
write.csv(Ftable,file = 'gene.list_ordered.csv')
