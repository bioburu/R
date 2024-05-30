library(org.Mm.eg.db)
li
retrieved <- AnnotationDbi::select(org.Mm.eg.db,
                                   keytype="GOALL",
                                   keys="GO:0030054",
                                   columns=c('SYMBOL','GENENAME','ENTREZID'))
