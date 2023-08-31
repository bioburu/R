reg<-FetchData(data,vars = c('ident',list),slot = 'counts')
ident<-as.data.frame(as.character(reg$ident))
reg<-reg[,-1]
reg<-t(reg)
reg<-as.data.frame(cpmNormalization(reg))
reg<-t(reg)
reg<-as.data.frame(cbind(ident,reg))
colnames(reg)[1]<-'ident'
