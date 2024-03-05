#BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
#BiocManager::install('biovizBase')
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
fragment<-file.path('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/geo_results/GSM6736815_IL1218_Neo_fragments.tsv.gz')
total_counts <- CountFragments(fragment)
head(total_counts)
summary(total_counts$frequency_count)
cutoff <- 1000 
barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB
head(barcodes)
frags <- CreateFragmentObject(path = fragment, cells = barcodes)
frags
peaks <- CallPeaks(frags,macs2.path = '/home/em_b/anaconda3/bin/macs3')
peaks
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
head(counts)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = frags,
  min.cells = 10,
  min.features = 200,
  genome = 'mm10')
data <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'ATAC',
  project = 'scATAC')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
annotations
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
annotations
Annotation(data) <- annotations
#---------------------------------------------------------------------------
data
data[['ATAC']]
granges(data)
head(data)
summary(data$nCount_ATAC)
VlnPlot(data,features=c('nCount_ATAC','nFeature_ATAC'),
        pt.size=0.1,
        ncol=2,
        cols = 'skyblue')
data<-subset(data,subset=nCount_ATAC>1000)
VlnPlot(data,features=c('nCount_ATAC','nFeature_ATAC'),
        pt.size=0.1,
        ncol=2,
        cols = 'skyblue')
cat('subset out all TSS scores below 4')
data <- NucleosomeSignal(object = data)
data <- TSSEnrichment(data)
DensityScatter(data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
head(data)
VlnPlot(object = data,
        features = c('nCount_ATAC',
                     'nFeature_ATAC',
                     'nucleosome_signal',
                     'nucleosome_percentile',
                     'TSS.enrichment','TSS.percentile'),
        pt.size = 0.1,
        ncol = 5,
        cols = 'skyblue')
data$high.tss <- ifelse(data$TSS.enrichment > 4, 'High', 'Low')
table(data$high.tss)
#TSSPlot(data, group.by = 'high.tss')# + NoLegend()
data<-subset(data,subset=TSS.enrichment>4)
VlnPlot(object = data,
        features = c('nCount_ATAC',
                     'nFeature_ATAC',
                     'nucleosome_signal',
                     'nucleosome_percentile',
                     'TSS.enrichment','TSS.percentile'),
        pt.size = 0.1,
        ncol = 5,
        cols = 'skyblue')
#-----------------------------------------------------------------------------------
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q0')
top500<-as.data.frame(FindTopFeatures(data[['ATAC']][])[1:500,])
View(top500)
data
data <- RunSVD(data) 
DepthCor(data)
cat('remove highly correlated dimensions')
data <- RunUMAP(object = data, reduction = 'lsi', dims = 3:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims =3:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3)
DimPlot(object = data, label = TRUE,pt.size = 2)
cat('cleaned ATACseq data is now ok for gene activity analysis')
gene.activities <- GeneActivity(data)
data[['correlated_gene_activity']] <- CreateAssayObject(counts = gene.activities)
cat('transfer to gene activity layer')
DefaultAssay(data)<-'correlated_gene_activity'
data
data <- NormalizeData(
  object = data,
  assay = 'correlated_gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_ATAC)
)
data<-FindVariableFeatures(data)
data<-ScaleData(data)
top1000 <- head(VariableFeatures(data), 1000)
top1000
FeaturePlot(
  object = data,
  features = c('Ptprc','Cd19','Cd22','Ncam1','Cd8a','Ifng','Il4','Cd3d','Il12rb1'),
  pt.size = 1,
  max.cutoff = 'q95',
  ncol = 3)
x<-FindMarkers(data, ident.1 = '1', ident.2 = '0', 
               features = c(top1000),test.use='negbinom')
VlnPlot(data, features = c('Ptprc','Cd19','Cd22','Ncam1','Cd8a','Ifng','Il4','Cd3d','Il12rb1'),cols = c())
VlnPlot(data, features = c(row.names(x)[1:12]),cols = c())
VlnPlot(data, features = c(row.names(x)[13:24]),cols = c())
VlnPlot(data, features = c(row.names(x)[25:36]),cols = c())
break
setwd('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/geo_results')
saveRDS(data, file = "geo.rds")
data<-readRDS('geo.rds')

