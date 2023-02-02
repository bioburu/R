library("GEOquery")
library("edgeR")
library("ggpubr")
library("ggthemes")
library("ggrepel")
library("RColorBrewer")
geo <-getGEO("GSE34860")
length(geo)
geo<- geo[[1]]
experimentData(geo)
lib <- pData(geo)
meta <- lib[,-c(1,3:40)]
annota <- fData(geo)
df <- data.frame(exprs(geo))
df<-na.omit(df)
design_factors<- meta$`npm status:ch1`
dge <- DGEList(counts=df,group=design_factors,remove.zeros = TRUE)
keep <- filterByExpr(dge, group = dge$samples$group, min.count = 100)
dge_edit <- dge[keep,]
dge_edit <-calcNormFactors(dge_edit)
multidimensional_scaling_plot <- plotMDS(dge_edit,method = "logFC",col = as.numeric(dge_edit$samples$group))
legend("topright",as.character(unique(dge_edit$samples$group)),cex=1,col=1:2,pch=20)
naive_exact_dispersions <- estimateDisp(dge_edit)
glm_dispersions <-estimateGLMCommonDisp(naive_exact_dispersions)
glm_dispersions <-estimateGLMTagwiseDisp(glm_dispersions) 
glm_dispersions <-estimateGLMTrendedDisp(glm_dispersions)
plotBCV(glm_dispersions) 
exact_test <- exactTest(naive_exact_dispersions, pair = c('wild-type','mutated'))
matrix <- model.matrix(~0+dge_edit$samples$group)
glmfit <-glmFit(glm_dispersions,matrix)
likelihood_ratio_test <- glmLRT(glmfit, contrast = c(1,-1))
quassi_likelihood_fit <- glmQLFit(glm_dispersions, matrix)
quassi_likelihood_Ftest <- glmQLFTest(quassi_likelihood_fit, contrast = c(1,-1))
topTags(exact_test, adjust.method = "BH", n = 10)
topTags(likelihood_ratio_test, adjust.method ="BH", n = 10)
topTags(quassi_likelihood_Ftest, adjust.method = "BH", n = 10)
DE_QLFtest <- decideTests(quassi_likelihood_Ftest, adjust.method = "BH", FDR = 0.05)
DE_tags <- row.names(dge_edit)[as.logical(DE_QLFtest)]
plotSmear(quassi_likelihood_Ftest,de.tags = DE_tags)
abline(h = c(-2,2), col = "blue")
plot(quassi_likelihood_Ftest$table$logFC, -log10(quassi_likelihood_Ftest$table$PValue),
     xlab = "logFC", ylab = "-log10(p.value)",
     cex = 0.5, pch = 5, col = "black")
points(quassi_likelihood_Ftest$table$logFC[quassi_likelihood_Ftest$table$PValue<0.000068],
       -log10(quassi_likelihood_Ftest$table$PValue[quassi_likelihood_Ftest$table$PValue<0.000068]),
              cex = 0.5, pch = 5, col = "red")
abline(v = 0, h = -log10(0.000068), lty = "dashed", col = "blue")
QLF_top50 <- data.frame(topTags(quassi_likelihood_Ftest, adjust.method = "BH", n = 50 ))
ID <- rownames(QLF_top50)
table <- cbind(ID, QLF_top50, stringsAsFactors = FALSE)
Ftable <- annota[annota$ID%in%c(table$ID),]
signif_gene_list <- merge(table, Ftable, by = "ID")
volcano_plot <- ggplot(signif_gene_list, aes(logFC, -log(FDR, 10)))+ 
  geom_point(pch = 5, size = 6, col = "grey")+
  geom_vline(xintercept = c(-1,1), col = "red")+
  geom_hline(yintercept = 5, col = "red")
volcano_plot <- volcano_plot + theme_few() + scale_colour_few()
gene_label <- signif_gene_list$`Gene Symbol`
volcano_plot <- volcano_plot +
  geom_label_repel(data = signif_gene_list,
                   mapping = aes(logFC, -log(FDR,10),label = gene_label),
                   size = 6, force = 2, nudge_x = -0.1, nudge_y = 0)+
  xlab(expression("log"[2]*"FC"))+
  ylab(expression("-log"[10]*"FDR"))
print(volcano_plot)
ID <-row.names(df)
df <- cbind(ID, df, stringsAsFactors = FALSE)
htmp_raw_counts <- df[df$ID%in%c(signif_gene_list$ID),]
gene_names <- annota[annota$ID%in%c(htmp_raw_counts$ID),]
gene_names <- gene_names[,-c(2:10,12:16)]
htmp_raw_counts <- merge(gene_names, htmp_raw_counts,by = "ID")
#break
#need to get rid of any duplicates if needed
#htmp_raw_counts[50,2] <- "MANIA1-2"
#rownames(htmp_raw_counts) <- htmp_raw_counts$`Gene Symbol`
#htmp_raw_counts <- htmp_raw_counts[,-c(1:2)]
#write.csv(htmp_raw_counts,file = "edit.csv")
#getwd()
htmp_raw_counts <- read.csv("edit.csv")
rownames(htmp_raw_counts) <- htmp_raw_counts$X
htmp_raw_counts <- htmp_raw_counts[,-1]           
htmp_raw_counts <- t(htmp_raw_counts)
htmp_cpm <- cpm(htmp_raw_counts)
htmp_logcpm <- cpm(htmp_raw_counts, log = TRUE)
color<-colorRampPalette(brewer.pal(8,"Blues"))(25)
heatmap(htmp_logcpm,scale='column',cexRow = 0.7, 
        cexCol = 0.7,col=color)
htmp_cpm <- cbind(meta, htmp_cpm, stringsAsFactors = FALSE)
htmp_cpm <- htmp_cpm[,-c(1,3)]
colnames(htmp_cpm)[1] <- "Genotype"
write.csv(htmp_cpm, file = "htmp_cpm.csv")
#htmp_cpm <- read.csv("htmp_cpm.csv")
#View(htmp_cpm)
dev.off()
device = "RStudioGD"
sessionInfo()
dev.cur()
