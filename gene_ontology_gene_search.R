pathway<-select(org.Mm.eg.db,
                keytype = 'GOALL',
                keys = 'GO:0007224',
                columns = c('SYMBOL','GENENAME','ENTREZID'))
list<-pathway$SYMBOL
list<-unique(list)

#----remove na's to construct seurat object 
gnps<-gnps[c(list),]
gnps<-rmdup(gnps)
gnps<-na.omit(gnps)
data <- CreateSeuratObject(counts=gnps,project = 'mb_gnp_hedgehog')
