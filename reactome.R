library(biomaRt)
library(celldex)
library(dplyr)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
geneid <- row.names(DEG)
head(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'human')
cat('Reactome pathway DEG annotations')
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways',
        label_format=50)

edox <- setReadable(reactome, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.5,
         circular=TRUE,
         cex_gene=0.1,
         showCategory = 10,
         node_label_size=1,
         cex_label_category=1)
