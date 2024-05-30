#-------------for Seurat object
FeatureScatter(data, feature1 = "Foxo1", feature2 = 'Dio2',cols = c('skyblue','red','grey','black'),pt.size = 5,
               shuffle = TRUE,seed = 123)
#---for data frame with variables as colnames (transposed)
library(ggpubr)
matrix<-read.csv('/home/em_b/work_stuff/chipseq/correlations/for_correlations.csv',
                 row.names = 1)
t<-data.frame(t(matrix))
colnames(t)
ggscatter(t,
          x='Neurod2',y='Adcy1',
          add='reg.line',
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = 'spearman',
          xlab='Neurod2',
          ylab = 'Adcy1',
          shape = 24,
          size=4,
          cor.coef.size = 4,
          ggtheme = theme_minimal(),
          label=row.names(t),
          font.label = c(10,'bold.italic','blue'),
          add.params = list(color='red',fill='lightgrey'))

#--------------for scRNAseq
df<-FetchData(data,vars = c('ident',top10000),layer = 'counts')
df<-df[,-1]
str(df)
cor<-data.frame(cor(df))
str(cor)
dio2<-subset(cor,GFAP>0.05)
heatmap(cor,scale = 'column')

write.csv(dio2,file = 'test_cor.csv')
