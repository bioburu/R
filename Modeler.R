x<-FindMarkers(data, ident.1 = 'LAC', ident.2 = 'adjacent', 
               features = c(top1000),logfc.threshold=1,min.pct1=1,
               max.pct2=0.0001,only.pos = TRUE)
write.csv(x,file = 'lung.pt1.csv')
VlnPlot(data, features = c(row.names(x)[1:3]),cols = c('grey','red'))
