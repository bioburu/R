x<-FindMarkers(data, ident.1 = 'HCC', ident.2 = 'adjacent', 
            features = c(gene_list),logfc.threshold=1,min.pct1=1,
            max.pct2=0.0001,only.pos = TRUE)
