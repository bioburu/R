library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
library(pheatmap)
exosomes<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/exosomes/csv_results/exosomes_biomart.csv',row.names = 1)
exosomes<-exosomes[,-c(1:8)]
exosomes<-round(exosomes)
#-------
tissues<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results/bulk_tissues_biomart.csv',row.names = 1)
tissues<-tissues[,-c(1:8)]
tissues<-round(tissues)
#------=
matrix<-cbind(tissues,exosomes)
matrix<-matrix[,-c(13:15)]
colnames(matrix)<-c('A.cereb.tiss_1','A.cereb.tiss_2','A.cereb.tiss_3',
                    'T.cereb.tiss_1','T.cereb.tiss_2','T.cereb.tiss_3',
                    'A.cereb.exo_1','A.cereb.exo_2','A.cereb.exo_3',
                    'A.pla.exo_1','A.pla.exo_2','A.pla.exo_3',
                    'T.pla.exo_1','T.pla.exo_2','T.pla.exo_3',
                    'T.cereb.exo_1','T.cereb.exo_2','T.cereb.exo_3')

data <- CreateSeuratObject(counts=matrix,project = 'exosomes')
gc()
table(data@active.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
top1000
gc()
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
gc()
data <- ScaleData(data, features = all.genes)
gc()
dim(data)
data <- RunPCA(data,npcs = 17, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:17, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
break 
VlnPlot(data,features = c('Atoh1','Eya1','Hhip','Pdlim3','Sfrp1','Mki67'),idents = c(),cols = c())
DimPlot(data,reduction = 'pca',pt.size = 4,cols = c('yellow','grey','green','red','skyblue','orange'),
        label = TRUE,repel = TRUE,label.box = TRUE,alpha = 3)
break 
VlnPlot(data,features = c('Rpl15.ps6','Rps2.5','Ncoa4.ps'),layer = 'counts',cols = c())
FeatureScatter(data,feature1 = 'Ncoa4.ps',feature2 = 'Rps2.5',
               shuffle = TRUE,seed = 123,pt.size = 7,slot = 'counts',
               plot.cor=FALSE,cols = c('yellow','grey','green','red','skyblue','orange'))
df<-FetchData(data,vars = c('ident','Rpl15.ps6','Rps2.5','Ncoa4.ps','Ezh2','Eya1'),layer = 'counts')
library(ggpubr)
ggscatter(df,x='Rpl15.ps6',y='Ncoa4.ps',
          add='reg.line',conf.int=TRUE,
          cor.coef = TRUE,cor.method = 'spearman',
          shape = 21,size = 6)
ggscatter(df,x='Rpl15.ps6',y='Rps2.5',
          add='reg.line',conf.int=TRUE,
          cor.coef = TRUE,cor.method = 'spearman',
          shape = 21,size = 6)
ggscatter(df,x='Ncoa4.ps',y='Rps2.5',
          add='reg.line',conf.int=TRUE,
          cor.coef = TRUE,cor.method = 'spearman',
          shape = 21,size = 6)
ggscatter(df,x='Ezh2',y='Eya1',
          add='reg.line',conf.int=TRUE,
          cor.coef = TRUE,cor.method = 'spearman',
          shape = 21,size = 6)
#-----------
break 
df<-FetchData(data,vars = c('ident','Rpl15.ps6','Rps2.5','Ncoa4.ps','Ezh2','Eya1'),layer = 'counts')

df$ident<-ifelse(df$ident=='T.pla.exo',1,0)
sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.5,0.5))
train  <- df[sample, ]
table(train$ident)
test   <- df[!sample, ]
table(test$ident)
model<-glm(train$ident~train$Rps2.5+train$Rpl15.ps6,family='binomial')
summary(model)
str(test)
predicted<-predict(model, test,type = 'response')
actual<-test$ident
predicted<-round(predicted)
predicted<-as.factor(predicted)
actual<-as.factor(actual)
confusion_matrix<-caret::confusionMatrix(predicted,actual,positive='1')
confusion_matrix
break 
#---------------
predicted<-predict(model, test,type = 'response')
predicted<-round(predicted)
actual<-test$ident

plot.roc(actual,predicted,percent=TRUE)
fourfoldplot(as.table(confusion_matrix))
break 
#-----A.cereb vs T.tiss
deseq_genes<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/exosomes/csv_results/exosomes_tissues_deseq2.csv')
deseq_genes<-deseq_genes$R_trackingID
deseq_genes
VlnPlot(data,features = c(deseq_genes[1:12]),idents = c('A.cereb','T.tiss'),cols = c('grey','red'))
VlnPlot(data,features = c(deseq_genes[13:24]),idents = c('A.cereb','T.tiss'),cols = c('grey','red'))
break 

#----A.pla vs t.pla
deseq_genes<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/exosomes/csv_results/exosomes_plasma_deseq2.csv')
deseq_genes<-deseq_genes$R_trackingID
deseq_genes
VlnPlot(data,features = c(deseq_genes[1:12]),idents = c('A.pla','T.pla'),cols = c('grey','red'))
VlnPlot(data,features = c(deseq_genes[13:24]),idents = c('A.pla','T.pla'),cols = c('grey','red'))

