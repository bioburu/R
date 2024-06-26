neuron_differentiation<-select(org.Mm.eg.db,
                               keytype = 'GOALL',
                               keys = 'GO:0030182',
                               columns = c('SYMBOL','GENENAME','ENTREZID'))
list<-neuron_differentiation$SYMBOL
