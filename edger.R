library(edgeR)
matrix<-read.csv('/home/bulk_tissues_biomart.csv',row.names = 1)
matrix<-matrix[,-c(1:8)]
colnames(matrix)<-c('Norm_1','Norm_2','Norm_3','MB_1','MB_2','MB_3')
condition <- factor(c(1,1,1,2,2,2))
dge <- DGEList(counts=matrix,group=condition)#,remove.zeros = TRUE)
keep<-filterByExpr(dge,group = condition)
DGE <- dge[keep,]
DGE <- normLibSizes(DGE)
DGE <- estimateDisp(DGE)
plotBCV(DGE)
fit <- glmQLFit(DGE)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)
df<-data.frame(qlf)
df<-df[order(df$PValue, decreasing=FALSE),]
df<-subset(df,PValue< 0.05)
#setwd('/home/csv_results')
#write.csv(df,file = 'bulk_tissues_edgeR.csv')
