library(sva)
matrix<-read.csv('/home/em_b/Desktop/chip_validations_final/in_house/rna_count_matrix.csv',
                 row.names = 1)
str(matrix)
matrix<-as.matrix(matrix)
batch<-c(1,1,1,2,2,2)
adjusted<-ComBat_seq(matrix,
                     batch = batch,
                     group = NULL)
adjusted
