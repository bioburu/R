df<-FetchData(data,vars = c('ident',top10000),layer = 'counts')
df<-df[,-1]
str(df)
cor<-data.frame(cor(df))
str(cor)
dio2<-subset(cor,GFAP>0.05)
heatmap(cor,scale = 'column')

write.csv(dio2,file = 'test_cor.csv')
