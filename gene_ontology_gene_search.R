pathway<-select(org.Mm.eg.db,
                keytype = 'GOALL',
                keys = 'GO:0030182',
                columns = c('SYMBOL','GENENAME','ENTREZID'))
list<-pathway$SYMBOL
