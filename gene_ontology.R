library(org.Hs.eg.db)
library(clusterProfiler)
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")# <--change to biological process, molecular funcion, or cellular component 
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        label_format=50)
