library(org.Mm.eg.db)

retrieved_genes <- AnnotationDbi::select(org.Mm.eg.db,
                                   keytype="GOALL",
                                   keys="GO:0030054",
                                   columns=c('SYMBOL','GENENAME','ENTREZID'))
list<-retrieved_genes$SYMBOL
list<-unique(list)
list
