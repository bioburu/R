x<-FindMarkers(data, ident.1 = '', ident.2 = '', 
               features = c(),logfc.threshold=1,min.pct1=1,
               max.pct2=0.0001,only.pos = TRUE)
write.csv(x,file = '')
VlnPlot(data, features = c(row.names(x)[1:3]),cols = c('grey','red'))
