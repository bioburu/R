library(pheatmap)
library(RColorBrewer)
matrix<-read.csv('/home/em_b/Desktop/chip_validations_final/in_house/rna_count_matrix.csv',
                 row.names = 1)
pheatmap(matrix[c('Neurod1','Insm1','Runx1t1','Tcf3'),],
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 17,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100))
#-------------------------------------------------------------------------------
library(rgl)
matrix<-read.csv('/home/em_b/Desktop/chip_validations_final/in_house/neurod_genes_3d.csv')
str(matrix)
pc <- princomp(matrix[,c(2,6,7,8)], cor=TRUE, scores=TRUE)
summary(pc)
plot(pc,type="lines")
biplot(pc)
setupKnitr(autoprint = TRUE)
plot3d(pc$scores,
       col=row.names(matrix),
       type='p',
       size = 10)

text3d(pc$scores,
       texts=matrix$condition,
       size=50)
text3d(pc$loadings, texts=rownames(pc$loadings), col="black")
coords <- NULL
for (i in 1:nrow(pc$loadings)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc$loadings[i,1:3]))
}
lines3d(coords, col="red", lwd=4)


https://planspacedotorg.wordpress.com/2013/02/03/pca-3d-visualization-and-clustering-in-r/
