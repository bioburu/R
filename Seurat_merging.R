library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/amp_prog/Downloads")
features_path <- "GSE130116_features.tsv.gz"
barcodes_path <- "GSE130116_barcodes.tsv.gz"
matrix_path <- "GSM3732336_ETV001_NYU_DIAGNOSIS.matrix.mtx.gz"
m1 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
obj1<-CreateSeuratObject(counts = m1,project = 'Diagnosis')
features_path <- "GSE130116_features.tsv.gz"
barcodes_path <- "GSE130116_barcodes.tsv.gz"
matrix_path <- "GSM3732337_ETV001_NYU_RELAPSE.matrix.mtx.gz"
m2 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
obj2<-CreateSeuratObject(counts = m2,project = 'Relapse')
data <- merge(obj1, y = obj2, add.cell.ids = c("Diagnosis", "Relapse"), project = "rrBALL")
break 
