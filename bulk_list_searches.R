library(dplyr)
setwd('/home/amp_prog/Downloads')
list1<-read.csv('milr1.cerv.1000.csv')
list2<-read.csv('milr1.ESCC.1000.csv')
list3<-read.csv('milr1.CRC.1000.csv')
joined<-semi_join(list1, list2, by = c('x'))
results<-semi_join(joined, list3, by = c('x'))
write.csv(results,file = 'results.csv')

###---to see significance of search results
setwd('/home/amp_prog/Downloads')
results<-read.csv('RESULTS.csv')
list<-results$X
setwd('/home/amp_prog/Downloads')
write.csv(x,file = 'milr1.CRC.1000.csv')
