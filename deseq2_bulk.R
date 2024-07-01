library(DESeq2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
setwd('/home/em_b/work_stuff/FCCC/T3_RNAseq/bulk_RNAseq_raw_files')
matrix<-read.csv('biomart.mm39.csv')
str(matrix)
#---if subsetting genes
list<-c('Neurod1','Cntn2','Map1a','Rbfox3',
        'Nefl','Fyn','Reln','Unc13a',
        'Ppp2r5b','Trim67','Cnr1','Kidins220',
        'Dpysl3','Rims1','Atp8a2','Mapt',
        'L1cam','Apbb1','Ss18l1',
        'Rufy3','Ndrg4','Pacsin1','Fkbp1b',
        'Neurod2','Sez6',
        'Map1b')
matrix<-subset(matrix, subset = Gene %in% list)
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-1]
head(matrix)
tester<-matrix[,-c(1:9,19:27)]
#tester<-tester[,-c(3,5,8)]
tester<-cbind(tester,tester[,c(4:6)])
tester<-tester[,-c(4:6)]
colnames(tester)<-c('GFP_overexp','GFP_overexp.1','GFP_overexp.2',
                    'shNeurod1','shNeurod1.1','shNeurod1.2',
                    'Neurod1_overexp','Neurod1_overexp.1','Neurod1_overexp.2')
#----deseq did not work
names<-colnames(tester)
names
condition<-c('A','A','A','A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
#-----------------------------------------------------------------------------
tester<-round(tester)
colnames(tester)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = tester,colData = coldata,design = ~ condition)
deseq
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
x<-data.frame(results)
x<-x[order(x$pvalue, decreasing=FALSE),]
x<-subset(x,padj< 0.05)
summary(x)
upreg<-subset(x,log2FoldChange> 0)
df<-cbind(row.names(upreg),upreg)
colnames(df)[1]<-'Gene'
#------------------------------------
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
dim(df)
geneid <- df$Gene
head(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('external_gene_name','entrezgene_id'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
#-------------------------------------------------------------------
genes$entrezgene_id
GO_result <- enrichGO(gene = genes$entrezgene_id,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
GO_result
df<-data.frame(GO_result)
break 
colnames(tester)<-c('GFP_overexp','GFP_overexp','GFP_overexp',
                    'shNeurod1','shNeurod1','shNeurod1',
                    'Neurod1_overexp','Neurod1_overexp','Neurod1_overexp')
pheatmap(tester[c('Neurod1','Cntn2','Map1a','Rbfox3',
                  'Nefl','Fyn','Reln','Unc13a',
                  'Ppp2r5b','Trim67','Cnr1','Kidins220',
                  'Dpysl3','Rims1','Atp8a2','Mapt',
                  'L1cam','Apbb1','Ss18l1',
                  'Rufy3','Ndrg4','L1cam','Pacsin1',
                  'Cxcl12','Fkbp1b','Neurod2','Sez6',
                  'Map1b'),],
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 17,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100))
