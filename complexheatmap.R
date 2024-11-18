library(DESeq2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
library(AnnotationDbi)
library(EnhancedVolcano)
library(stringr)
library(fgsea)
library(Seurat)
library(ComplexHeatmap)
setwd('/home/')
matrix<-read.csv('GSE224974_mm39_raw_count_matrix.csv',
                 row.names = 1)
head(matrix)
list<-row.names(matrix)
list<-list[!grepl('Gm',list)]
list<-list[!grepl('Rik',list)]
matrix<-subset(matrix,row.names(matrix)%in%list)
matrix<-matrix[,-c(7:9)]
matrix<-round(matrix)
colnames(matrix)<-c('Zhr_1','Zhr_2','Zhr_3',
                    'T3_1','T3_2','T3_3')
data <- CreateSeuratObject(counts=matrix,project = 'mm39')
data@active.ident
table(data@active.ident)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
all.genes <- rownames(data)
data <- ScaleData(data,
                  do.center = FALSE,
                  features = all.genes)
data <- RunPCA(data,npcs = 5, features = VariableFeatures(object = data))
VlnPlot(data,
        features = c('Neurod1'),
        pt.size = 0.1,
        cols = c('red','skyblue'),
        layer='counts',
        raster = TRUE)
#----smoothened signaling pathway ----------------------------------------------
pathway<-select(org.Mm.eg.db,
                keytype = 'GOALL',
                keys = 'GO:0007224',
                columns = c('SYMBOL','GENENAME','ENTREZID'))
list<-pathway$SYMBOL
matrix_subset<-subset(matrix, subset = row.names(matrix) %in% list)
names<-colnames(matrix_subset)
names
condition<-c('A','A','A',
             'B','B','B')
type<-c('paired','paired','paired',
        'paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
coldata
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
colnames(matrix_subset)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = matrix_subset,
                              colData = coldata,
                              design = ~ condition)
deseq
DE<-DESeq(deseq)
results<-results(DE)
results
EnhancedVolcano(results,
                lab = row.names(results),
                x='log2FoldChange',
                y='pvalue',
                title = 'GO:0007224. Smoothened signaling pathway',
                subtitle = '0hr vs T3 treatments',
                legendLabels = NULL,
                legendIconSize = -1,
                legendPosition = 'bottom',
                pCutoff = 0.05,
                FCcutoff = 0.1,
                shape = 1,
                ylim = c(0,14),
                xlim = c(-5,2.5))
results_df<-data.frame(results)
#-------------------------------------------------------------------------------
upreg<-subset(results_df,log2FoldChange> 0.1)
upreg<-subset(upreg,pvalue< 0.05)
pheatmap(data.matrix(matrix[c(row.names(upreg)),]),
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'green',
         fontsize = 10,
         legend = FALSE,
         display_numbers = TRUE,
         number_color = 'blue',
         fontsize_number = 15,
         main = 'GO:0007224:_Smoothened_signaling_pathway',
         color = colorRampPalette(rev(brewer.pal(n=15,name='RdYlBu')))(100))
test<-matrix[c(row.names(upreg)),]
test<-data.matrix(test, rownames.force = NA)
col_fun = circlize::colorRamp2(c(0, 100, 400, 800), c("white", "skyblue", "orange", "red"))
ht_opt$TITLE_PADDING = unit(c(10, 1), "mm")
Heatmap3D(test,
          bar_rel_width = 0.4,
          bar_rel_height = 0.7,
          bar_angle = 50,
          bar_max_length = unit(1,'cm'),
          name='Counts',
          column_title='GO:0007224:_Smoothened_signaling_pathway',
          col=col_fun,
          show_row_dend = FALSE,
          cluster_columns=FALSE,
          row_names_side = 'right',
          row_names_gp=gpar(fontsize=15),
          heatmap_legend_param = list(at = c(0, 400, 700, 1000), 
                                      labels = c('0','400','700','1000')))
#-------------------------------------------------------------------------------
downreg<-subset(results_df,log2FoldChange< -0.1)
downreg<-subset(downreg,pvalue< 0.05)
pheatmap(data.matrix(matrix[c(row.names(downreg)),]),
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'green',
         fontsize = 10,
         legend = FALSE,
         display_numbers = TRUE,
         number_color = 'blue',
         fontsize_number = 9,
         main = 'GO:0007224:_Smoothened_signaling_pathway',
         color = colorRampPalette(rev(brewer.pal(n=15,name='RdYlBu')))(100))
test<-matrix[c(row.names(downreg)),]
test<-data.matrix(test, rownames.force = NA)
col_fun = circlize::colorRamp2(c(0, 10, 100, 500), c("white", "skyblue", "orange", "red"))
ht_opt$TITLE_PADDING = unit(c(10, 1), "mm")
Heatmap3D(test,
          bar_rel_width = 0.4,
          bar_rel_height = 0.7,
          bar_angle = 50,
          bar_max_length = unit(1,'cm'),
          name='Counts',
          column_title='GO:0007224:_Smoothened_signaling_pathway',
          col=col_fun,
          show_row_dend = FALSE,
          cluster_columns=FALSE,
          row_names_side = 'right',
          row_names_gp=gpar(fontsize=14),
          heatmap_legend_param = list(at = c(0, 10, 100, 500), 
                                      labels = c('0','10','100','500'))) 
#----chromatin lock complex ----------------------------------------------
pathway<-select(org.Mm.eg.db,
                keytype = 'GOALL',
                keys = 'GO:0061793',
                columns = c('SYMBOL','GENENAME','ENTREZID'))
list<-pathway$SYMBOL
matrix_subset<-subset(matrix, subset = row.names(matrix) %in% list)
pheatmap(matrix_subset,
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'green',
         fontsize = 13,
         legend = FALSE,
         display_numbers = TRUE,
         number_color = 'blue',
         fontsize_number = 15,
         main = 'GO:0061793:_Chromatin_lock_complex',
         color = colorRampPalette(rev(brewer.pal(n=15,name='RdYlBu')))(100))
test<-data.matrix(matrix_subset, rownames.force = NA)
col_fun = circlize::colorRamp2(c(0, 10, 1800, 2500), c("white", "skyblue", "orange", "red"))
ht_opt$TITLE_PADDING = unit(c(10, 1), "mm")
Heatmap3D(test,
          bar_rel_width = 0.4,
          bar_rel_height = 0.7,
          bar_angle = 50,
          bar_max_length = unit(1,'cm'),
          name='Counts',
          column_title='GO:0061793:_Chromatin_lock_complex',
          col=col_fun,
          show_row_dend = FALSE,
          cluster_columns=FALSE,
          row_names_side = 'right',
          row_names_gp=gpar(fontsize=15),
          heatmap_legend_param = list(at = c(0, 10, 1800, 2500), 
                                      labels = c('0','10','1800','2500')))
#----chromatin silencing complex  ----------------------------------------------
pathway<-select(org.Mm.eg.db,
                keytype = 'GOALL',
                keys = 'GO:0005677',
                columns = c('SYMBOL','GENENAME','ENTREZID'))
list<-pathway$SYMBOL
matrix_subset<-subset(matrix, subset = row.names(matrix) %in% list)
pheatmap(matrix_subset,
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'green',
         fontsize = 12,
         legend = FALSE,
         display_numbers = TRUE,
         number_color = 'blue',
         fontsize_number = 15,
         main = 'GO:0005677:_Chromatin_silencing_complex',
         color = colorRampPalette(rev(brewer.pal(n=15,name='RdYlBu')))(100))
test<-data.matrix(matrix_subset, rownames.force = NA)
col_fun = circlize::colorRamp2(c(0, 5, 1500, 2000), c("white", "skyblue", "orange", "red"))
ht_opt$TITLE_PADDING = unit(c(10, 1), "mm")
Heatmap3D(test,
          bar_rel_width = 0.4,
          bar_rel_height = 0.7,
          bar_angle = 50,
          bar_max_length = unit(1,'cm'),
          name='Counts',
          column_title='GO:0005677:_Chromatin_silencing_complex',
          col=col_fun,
          show_row_dend = FALSE,
          cluster_columns=FALSE,
          row_names_side = 'right',
          row_names_gp=gpar(fontsize=15),
          heatmap_legend_param = list(at = c(0, 10, 1800, 2500), 
                                      labels = c('0','10','1800','2500')))
#-------------------------------------------------------------------------------
