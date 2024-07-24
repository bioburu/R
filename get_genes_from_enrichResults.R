reactome <- enrichPathway(upreg_matrix$entrezgene_id,
                          organism = 'mouse')
reactome_df<-data.frame(reactome)
reactome_df<-reactome_df[order(reactome_df$Count, decreasing=TRUE),]
pathways<-reactome_df$ID
reactome@result<-reactome@result[reactome@result$ID%in%pathways,]
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=15,
        label_format=50,
        title='Tu-AST and norm-AST Reactome pathways upregulated genes',
        showCategory=24)
View(data.frame(reactome))
target_genes<-reactome@result$geneID[9]
target_genes
target_genes<-unlist(strsplit(target_genes,split = '/'))
