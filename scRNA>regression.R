setwd('/home/amp_prog/Desktop/hope_model/GSE208653_hpv_cancer')
features_path <- 'GSM6360686_SCC_4.features.tsv.gz'
barcodes_path <- 'GSM6360686_SCC_4.barcodes.tsv.gz'
matrix_path <- 'GSM6360686_SCC_4.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
df<-as.data.frame(matrix[c('CD14','CD274'),])
df<-df[,c(1:500)]
df<-as.data.frame(t(df))
write.csv(df,file = 'cancer.csv')
