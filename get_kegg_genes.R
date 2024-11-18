matrix<-read.csv('/home/Desktop/prion_disease_kegg.csv',
                 header = FALSE)
v<-matrix$V1
library(stringr)
v
v2<-str_extract(v, "[^;]+")
v2<-na.omit(v2)
v2<-v2[-c(70,275:277)]
v2
cat(v2)
write.csv(v2,file = '/home/Desktop/prion_disease_kegg.csv')
