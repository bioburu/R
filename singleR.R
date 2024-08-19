library(celldex)
library(SingleR)
#--auto annotate cell clusters
surveyReferences()
searchReferences('mouse')
ref <- fetchReference("mouse_rnaseq", "2024-02-26")
ref
#--- convert seurat object into single cell experiment object 
rna.sce <- as.SingleCellExperiment(DietSeurat(rna))
rna.sce
ref.main <- SingleR(test = rna.sce,
                       assay.type.test = 1,
                       ref = ref,
                       labels = ref$label.main)
ref.fine <- SingleR(test = rna.sce,
                    assay.type.test = 1,
                    ref = ref,
                    labels = ref$label.fine)
table(ref.main$pruned.labels)
table(ref.fine$pruned.labels)
rna@meta.data$ref.main <- ref.main$pruned.labels
rna@meta.data$ref.fine <- ref.fine$pruned.labels
head(rna@meta.data)
rna <- SetIdent(rna, value = "ref.fine")
