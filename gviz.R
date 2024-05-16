library(biomaRt)
library(dplyr)
library(Gviz)
library(IRanges)
library(IRanges)
#setwd('/home/em_b/Downloads')
#matrix <- read.csv('matrix.csv')
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
filters
fm=Gviz:::.getBMFeatureMap()
fm['symbol']='ensembl_gene_id'
bm = BiomartGeneRegionTrack(chromosome='chr2', genome="mm39", 
                            start=79162835, end = 79382744, 
                            biomart=ensembl,
                            #filter=list("with_ox_refseq_mrna"=TRUE), 
                            size=4, name="RefSeq", utr5="red3", utr3="red3", 
                            protein_coding="black", col.line=NULL, cex=7,
                            collapseTranscripts="longest",
                            featureMap=fm)
bm
AT=GenomeAxisTrack()
plotTracks(c( bm, AT),
           from=79162835, to=79382744,
           transcriptAnnotation="symbol", window="auto", 
           cex.title=1, fontsize=10 )

#----------------
