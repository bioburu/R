library(plotly)
#all.genes <- rownames(data)
#data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
#DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
#ElbowPlot(data)
#data<-JoinLayers(data)
data <- RunUMAP(data,dims = 1:15,n.components = 3L)
head(data@meta.data)
plot3d1 <- FetchData(data, vars = c("umap_1", "umap_2", "umap_3", "orig.ident","cell_types"))
plot3d1$label <- paste(plot3d1$orig.ident)
fig <- plot_ly(data = plot3d1, 
               x = ~umap_1, y = ~umap_2, z = ~umap_3, 
               color = ~orig.ident, 
               colors = c('viridis'),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 2, width=2),
               text=~label,
               hoverinfo="text")%>%layout(title='GSE217511-UMAP_of_17wks_and_20wks')
fig
