#------Remove duplicate genes by p_value
matrix<-matrix%>%
  group_by(Gene)%>%
  arrange(p_value)%>%
  slice(1)
