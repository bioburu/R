library(dplyr)
setwd('/home/amp_prog/Downloads')
list1<-read.csv('cervical.csv')
list2<-read.csv('CRC.csv')
list3<-read.csv('ESCC.csv')
joined<-semi_join(list1, list2, by = c('x'))
results<-semi_join(joined,list3,by=c('x'))
write.csv(results,file = 'results.csv')

###---to see significance of search results
setwd('/home/amp_prog/Downloads')
results<-read.csv('results.csv')
list<-results$x
x<-FindMarkers(data, ident.1 = 'ESCC', ident.2 = 'adjacent', features = c(list))
write.csv(x,file = 'sig.ESCC.csv')
