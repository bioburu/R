library(volcano3D)
library(knitr)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
library(AnnotationDbi)
library(EnhancedVolcano)
library(stringr)
library(ComplexHeatmap)
setwd('/home/em_b/work_stuff/*Differentiation_study')
matrix<-read.csv('mm39_raw_count_matrix.csv',
                 row.names = 1)
head(matrix)
list<-row.names(matrix)
list<-list[!grepl('Gm',list)]
list<-list[!grepl('Rik',list)]
matrix<-subset(matrix,row.names(matrix)%in%list)
matrix<-matrix[,c(1:3,7:12)]
matrix<-round(matrix)
colnames(matrix)<-c('Zhr_1','Zhr_2','Zhr_3',
                    'PBS_1','PBS_2','PBS_3',
                    'T3_1','T3_2','T3_3')
names<-colnames(matrix)
names
condition<-c('0hr','0hr','0hr',
             'PBS','PBS','PBS',
             'T3','T3','T3')
type<-c('paired','paired','paired',
        'paired','paired','paired',
        'paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
coldata
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
colnames(matrix)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = matrix,
                              colData = coldata,
                              design = ~ condition)
deseq
DE<-DESeq(deseq)
DE
LRT<-DESeq(deseq,
          test = c('LRT'),
          reduced = ~ 1,
          parallel = TRUE)
res<-deseq_polar(DE,LRT,'condition')
str(res)
volcano3D(
  res,
  type = 2,
  label_rows = c('Neurod1','Gli2','Klf9'),
  label_size = 14,
  arrow_length = 120,
  colour_code_labels = TRUE,
  label_colour = "black",
  grid_colour = "grey80",
  grid_width = 2,
  grid_options = NULL,
  axis_colour = "red",
  axis_width = 2,
  marker_size = 5,
  marker_outline_width = 1,
  marker_outline_colour = "grey",
  z_axis_title_offset = 1.2,
  z_axis_title_size = 12,
  z_axis_angle = 0.5,
  radial_axis_title_size = 20,
  radial_axis_title_offset = 1.2,
  xy_aspectratio = 3,
  z_aspectratio = 1,
  camera_eye = list(x = 0.9, y = 0.9, z = 0.9))
radial_ggplot(res,
              marker_size = 2.3,
              legend_size = 10) +
  theme(legend.position = "right")
