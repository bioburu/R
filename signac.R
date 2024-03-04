#BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
#BiocManager::install('biovizBase')
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
counts<-readRDS('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/SRR22324456/batch1/sparse_matrix.rds')
colnames(counts)
row.names(counts)
cat('sort -k1,1 -k2,2n -k3,3n fragments.bed >fragments.sorted.bed')
cat('bgzip fragments.sorted.bed')
cat('tabix -p bed fragments.sorted.bed.gz')
cat('convert all bed files to tsv')
fragment<-file.path('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/SRR22324456/batch1/fragments.sorted.tsv.gz')
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = fragment,
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
break
saveRDS(data, file = "batch1.rds")
data<-readRDS('batch1.rds')
