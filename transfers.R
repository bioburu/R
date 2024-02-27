head(rna)
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna)
DimPlot(rna, reduction = "umap",dims = c(1,2))+ggtitle('scRNAseq')
cluster.markers <- FindMarkers(rna, ident.1 = 0)
head(cluster.markers,n=25)
break 
#---convert into Assay5 obj
#pbmc.rna[['RNA']]<-as(pbmc.rna[['RNA']],Class = 'Assay5')
# get most variable gene activity from the RNAseq file of the ATACseq file
gene.activities <- GeneActivity(atac, features = VariableFeatures(rna))

break 
