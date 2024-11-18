library(DESeq2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
library(EnhancedVolcano)
library(ReactomePA)
setwd('/home/AST_OPC_files')
matrix<-read.csv('matrix_mm39.csv')
str(matrix)
row.names(matrix)<-make.names(matrix$external_gene_name,unique = TRUE)
matrix<-matrix[,-c(1:2)]
head(matrix)
tester<-round(matrix)
astrocytes<-tester[,-c(5:8,11:12)]
names<-colnames(astrocytes)
names
condition<-c('A','A','A','A','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
#-----------------------------------------------------------------------------
colnames(astrocytes)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = astrocytes,colData = coldata,design = ~ condition)
deseq
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
EnhancedVolcano(results,
                lab = row.names(results),
                x='log2FoldChange',
                y='pvalue',
                title = 'Differential gene expressions of Tu-AST and norm-AST',
                subtitle = '',
                legendLabels = NULL,
                legendIconSize = -1,
                legendPosition = 'bottom',
                pCutoff = 0.01,
                FCcutoff = 1.5,
                shape = 1,
                ylim = c(0,100),
                xlim = c(-9,11))
#-----------------------------------------------------------------------
x<-data.frame(results)
x<-x[order(x$pvalue, decreasing=FALSE),]
x<-subset(x,pvalue< 0.01)
summary(x)
upreg<-subset(x,log2FoldChange> 1.5)
#------------------------------------------------------------------
matrix<-read.csv('matrix_mm39.csv')
upreg_matrix<-subset(matrix,matrix$external_gene_name%in%row.names(upreg))
reactome <- enrichPathway(upreg_matrix$entrezgene_id,
                          organism = 'mouse')
reactome_df<-data.frame(reactome)
reactome_df<-reactome_df[order(reactome_df$Count, decreasing=TRUE),]
pathways<-reactome_df$ID
reactome@result<-reactome@result[reactome@result$ID%in%pathways,]
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=15,
        label_format=50,
        title='Tu-AST and norm-AST Reactome pathways upregulated genes',
        showCategory=24)
#-------------------------------------------------------------------------------
go <- enrichGO(gene = upreg_matrix$entrezgene_id,
               OrgDb = org.Mm.eg.db,
               ont = "BP")
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=15,
        label_format=50,
        title='Tu-AST and norm-AST GO:Biological Process pathways upregulated genes',
        showCategory=30)
#--------------------------------------------------------------------
go <- enrichGO(gene = upreg_matrix$entrezgene_id,
               OrgDb = org.Mm.eg.db,
               ont = "CC")
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=15,
        label_format=50,
        title='Tu-AST and norm-AST GO:Cellular Component pathways upregulated genes',
        showCategory=30)
#-------------------------------------------------------------------------------
go <- enrichGO(gene = upreg_matrix$entrezgene_id,
               OrgDb = org.Mm.eg.db,
               ont = "MF")
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=15,
        label_format=50,
        title='Tu-AST and norm-AST GO:Molecular Function pathways upregulated genes',
        showCategory=30)

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
