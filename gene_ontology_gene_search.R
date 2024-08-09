pathway<-select(org.Mm.eg.db,
                keytype = 'GOALL',
                keys = 'GO:0007224',
                columns = c('SYMBOL','GENENAME','ENTREZID'))
list<-pathway$SYMBOL
list<-unique(list)
