library(biomaRt)
library(dplyr)
matrix <- read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results/bulk-tissues_quant.csv')
matrix<-matrix[,-c(2:3)]
#---set host to https://useast.ensembl.org or https://uswest.ensembl.org
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
dim(matrix)
geneid <- matrix$Gene
head(geneid)
dim(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_transcript_id_version','entrezgene_id','external_gene_name','chromosome_name','p_value','start_position','end_position','description'),
              filters = 'ensembl_transcript_id_version',
              values = geneid,
              mart = ensembl)
head(genes)
df <- matrix[matrix$Gene %in% unique(genes$ensembl_transcript_id_version),]
dim(df)
dim(matrix)
colnames(genes)[1]<- 'Gene'
Ftable <- merge(genes, df, by="Gene")
colnames(Ftable)[1] <- c('transcript_version')
colnames(Ftable)[3] <- c('Gene')
Ftable<-Ftable[!duplicated(Ftable$transcript_version),]
Ftable<-data.frame(Ftable)
row.names(Ftable)<-make.names(Ftable$Gene,unique = TRUE)
str(Ftable)
Ftable<-cbind(row.names(Ftable),Ftable)
colnames(Ftable)[4]<-'external_gene_name'
colnames(Ftable)[1]<-'R_trackingID'
break
setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results')
#write.csv(Ftable,file = 'bulk_tissues_biomart.csv')
library(DESeq2)
library(pheatmap)
library(dplyr)
matrix<-Ftable[,-c(1:9)]
colnames(matrix)<-c('Norm_1','Norm_2','Norm_3','MB_1','MB_2','MB_3')
matrix<-round(matrix)
head(matrix)
names<-colnames(matrix)
condition<-c('A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-coldata$names
coldata<-coldata[,-1]
head(coldata)
#-----------------------------------------------------------------------------
deseq<-DESeqDataSetFromMatrix(countData = matrix,colData = coldata,design = ~ condition)
deseq
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-2,2))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
x<-data.frame(results)
summary(x)
x<-na.omit(x)
#----only upregulated genes 
downreg<-subset(x,log2FoldChange< -1.5)
upreg<-subset(x,log2FoldChange> 1.5)
df<-rbind(upreg,downreg)

df<-cbind(row.names(df),df)
df$`row.names(df)`<-sub('\\..*','',df$`row.names(df)`)
colnames(df)[1]<-'external_gene_name'
df<-cbind(row.names(df),df)
colnames(df)[1]<-'R_trackingID'
biomart_output<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results/bulk_tissues_biomart.csv',row.names = 1)
df <- merge(df, biomart_output,by="R_trackingID")
df<-df[,-11]
colnames(df)[2]<-'external_gene_name'
df<-df[order(df$padj, decreasing=FALSE),]
df<-subset(df,padj< 0.05)
df<-df[,-12]
setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results')
#write.csv(df,file = 'bulk_tissues_deseq2.csv')
#--------------
library(pathview)
foldchanges = df$log2FoldChange
names(foldchanges) = df$entrezgene_id
head(foldchanges)
pv.out <- pathview(gene.data = foldchanges, pathway.id = "05022",
                   species = "mmu", out.suffix = "mmu_test")
break 
