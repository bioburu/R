#-------------for Seurat object
FeatureScatter(data, feature1 = "Foxo1", feature2 = 'Dio2',cols = c('skyblue','red','grey','black'),pt.size = 5,
               shuffle = TRUE,seed = 123)
#---for data frame with variables as colnames (transposed)
ggscatter(df,x='Foxo1',y='Dio2',
          add='reg.line',conf.int = TRUE,
          cor.coef = TRUE,cor.method = 'spearman',
          xlab='Foxo1',ylab = 'Dio2',shape = 21,size=2)

#--------------for scRNAseq
df<-FetchData(data,vars = c('ident',top10000),layer = 'counts')
df<-df[,-1]
str(df)
cor<-data.frame(cor(df))
str(cor)
dio2<-subset(cor,GFAP>0.05)
heatmap(cor,scale = 'column')

write.csv(dio2,file = 'test_cor.csv')
