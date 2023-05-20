setwd('')
features_path <- 'GSM5900215_1T_features.tsv.gz'
barcodes_path <- 'GSM5900215_1T_barcodes.tsv.gz'
matrix_path <- 'GSM5900215_1T_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
norm<-as.data.frame(matrix[c('CD274','EBNA1BP2'),])
factor<-rep(c(1),times=nrow(norm))
norm<-rbind(factor,norm)
row.names(norm)[1]<-'group'
#---------------------------------------------------------------------------
features_path <- 'GSM5900216_1N_features.tsv.gz'
barcodes_path <- 'GSM5900216_1N_barcodes.tsv.gz'
matrix_path <- 'GSM5900216_1N_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
norm2<-as.data.frame(matrix[c('CD274','EBNA1BP2'),])
factor<-rep(c(0),times=nrow(norm2))
norm2<-rbind(factor,norm2)
row.names(norm2)[1]<-'group'
out<-cbind(norm,norm2)
out<-t(out)
write.csv(out,file = 'out.csv')
