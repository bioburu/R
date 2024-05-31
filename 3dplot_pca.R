library(rgl)
matrix<-read.csv('test.csv')

pc <- princomp(matrix[,1:4], cor=TRUE, scores=TRUE)
summary(pc)
plot(pc,type="lines")
biplot(pc)
setupKnitr(autoprint = TRUE)
plot3d(pc$scores,
       col=matrix$sample_id,
       type='p',
       size = 10)

text3d(pc$scores,
       texts=rownames(matrix),
       size=50)
text3d(pc$loadings, texts=rownames(pc$loadings), col="black")
coords <- NULL
for (i in 1:nrow(pc$loadings)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc$loadings[i,1:3]))
}
lines3d(coords, col="grey", lwd=4)
