x<-FindMarkers(data, ident.1 = '', ident.2 = '', 
               features = c(),logfc.threshold=1.5,only.pos = TRUE,test.use='DESeq2')
VlnPlot(data, features = c(row.names(x)[1:12]),cols = c('grey','red'))
write.csv(x,file = '')
