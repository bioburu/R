both_common
both_common<-both_common[!grepl('Gm',both_common)]
both_common
both_common<-both_common[!grepl('Rik',both_common)]

#---for row.names
markers<-markers[-grep("Rik", row.names(markers)),]
markers<-markers[-grep("mt-", row.names(markers)),]
