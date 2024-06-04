library(ggvenn)
library(gplots)

x1<-read.csv('THRA3rep_0hr.csv')
x2<-read.csv('THRA3rep_B22.csv')
x3<-read.csv('THRA3rep_T3.csv')
x <- list(
  THRA_0hr = x1$SYMBOL, 
  THRA_B22 = x2$SYMBOL, 
  THRA_T3 = x3$SYMBOL)
x
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
print(v.table)
