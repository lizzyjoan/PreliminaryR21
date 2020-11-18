
---
  title: "Test"
author: John Doe
date: 22/Oct/2016
output: html_document
---

```{r}
library(projectR)
#for geo stuff
library(Biobase)
library(GEOquery)
library(limma)
#deseq for normalizing
BiocManager::install("DESeq2")
#gplots for NMF mpas
library(gplots)

#read in mouse data
mouse_data <- read.delim("fixed_raw_counts.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
mouse_data <- mouse_data[, -c(1:3)]
#mouse_list <- list(Type = c("WT_NC", "WT_NC", "WT_NC", "WT_HFD", "WT_HFD", "WT_HFD", "WT_HFD_AQ", "WT_HFD_AQ", "WT_HFD_AQ", "WT_NC_AQ", "WT_NC_AQ", "WT_NC_AQ", "HFD_WT", "HFD_WT", "HFD_WT", "HFD_NURR", "HFD_NURR", "HFD_NURR"))
#mouse_diet <- list(Diet = c("Normal", "Normal", "Normal", "High_Fat", "High_Fat", "High_Fat", "High_Fat", "High_Fat", "High_Fat", "Normal", "Normal", "Normal", "High_Fat", "High_Fat", "High_Fat", "High_Fat", "High_Fat", "High_Fat"))
#mouse_genotype <- list(Genotype = c("WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "NURR", "NURR", "NURR"))
#mouse_meta <- cbind(mouse_list, mouse_diet, mouse_genotype)
#mouse_meta <- as.data.frame(mouse_meta)
#mouse_meta <- unlist(mouse_meta)
mouse_info <- read.csv("mouse_info.csv", header = TRUE)

####DESeq2 vst for normalizing
#MOUSE
t_mouse_data <- t(mouse_info)
mouse_matrix <- cbind(t_mouse_data, mouse_info)
mouse_dds <- DESeqDataSetFromMatrix(mouse_data, colData = mouse_info, design = ~ AQ)
mouse_dds <- DESeq(mouse_dds)
mouse_res <- results(mouse_dds)
mouse_vsd <- vst(mouse_dds, blind = FALSE)
mouse_deseq2_pca <- plotPCA(mouse_vsd, intgroup = c("AQ", "Diet"))
mouse_deseq2_pca

#Human data read in
#human_metadata <- read.csv("human_metadata.csv", header = TRUE, sep = ",")
#human_metadata <- t(human_metadata)
human_data <- read.delim("human_fixed_raw_counts.txt", header = TRUE, row.names = 1, sep = "\t")
human_data <- human_data[, -c(1:3)]
#human_metadata <- read.table("Human_GSE108643_Cuffdiff_gene_differential_expression_testing.tabular.txt", fill = TRUE)
human_info <- read.csv("human_info.csv", sep = ",")
human_info <- human_info[c(1:58),]

#HUMAN normalization
#t_human_data <- t(human_info)
#human_matrix <- cbind(t_human_data, human_info)
human_dds <- DESeqDataSetFromMatrix(human_data, colData = human_info, design = ~ PrePost)
human_dds <- DESeq(human_dds)
human_res <- results(human_dds)
human_vsd <- vst(human_dds, blind = FALSE)
human_deseq2_pca <- plotPCA(human_vsd, intgroup = c("PrePost", "BMI"))
huma_deseq2_pca


####PCA calculation
#MOUSE
pc.mouse_data <- prcomp(assay(mouse_vsd))
pcVAR_mouse <- round(((pc.mouse_data$sdev)^2/sum(pc.mouse_data$sdev^2))*100,2)
dfPCA_mouse <- data.frame(pc.mouse_data$x)

#HUMAN
pc.human_data <- prcomp(assay(human_vsd))

#Obtaining PC's to project
#plot PCA
library(ggplot2)
setCOL <- scale_colour_manual(values = c("blue", "black", "red"), name = "Diet")
setFILL <- scale_fill_manual(values = c("blue", "black", "red"), guide = FALSE)
setPCH <- scale_shape_manual(values = c(23, 22, 25, 25, 21, 24), name = "Genotype")

pPCA <- ggplot(dfPCA_mouse, aes(x = PC1, y = PC2, colour = Diet, shape = Genotype, 
                          fill = Genotype)) + geom_point(aes(size = 1), alpha = .6) +
  setCOL + setPCH + setFILL +
  scale_size_area(breaks = c(2, 4, 6) ) +
  theme(legend.position = c(0.6,0.6), legend.justification = c(0,0),
        legend.direction = "horizontal",
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(vjust = 0, hjust = 0, face = "bold")) +
  labs(title = "PCA of Mouse RNAseq",
       x = paste("PC1 (", pcVAR[1], "% of variance)", sep = ""),
       y = paste("PC2 (", pcVAR[2], "% of variance)", sep = ""))
pPCA

#trying a simpler pca plot from https://www.datacamp.com/community/tutorials/pca-analysis-r?utm_source=adwords_ppc&utm_campaignid=1565261270&utm_adgroupid=67750485268&utm_device=c&utm_keyword=&utm_matchtype=b&utm_network=g&utm_adpostion=&utm_creative=295208661496&utm_targetid=aud-299261629574:dsa-429603003980&utm_loc_interest_ms=&utm_loc_physical_ms=9012546&gclid=CjwKCAiA4o79BRBvEiwAjteoYPiR4NccRCpN9162HbTz55R_J6Iq0GC_jjA7czwSoKSBUYqF8rM8NBoCO1cQAvD_BwE
ggbiplot(pc.mouse_data)

####
#trying different normalization 
library(Biobase)
library(GenomicRanges)
projection <- projectR(data = human_data, loadings = pc.mouse_data, full = TRUE)

human_matrix <- as.matrix(human_data)
human_se <- SummarizedExperiment(human_matrix)
human_se_2 <- as.
human_vsd <- assay(human_vsd)
AP.human_matrix <- CoGAPS(human_vsd, nIterations = 10)

AP.human_mean <- AP.human_matrix$Amean
pNMF <-heatmap.2(as.matrix(AP.human_matrix),col=bluered, trace='none',
                 distfun=function(c) as.dist(1-cor(t(c))) ,
                 cexCol=1,cexRow=.5,scale = "row",
                 hclustfun=function(x) hclust(x, method="average"))
plot(AP.human)
pNMF

saveRDS(AP.mouse, "AP.mouse")
saveRDS(AP.human, "AP.human")
saveRDS(human_svd <- human_svd, "human_svd")

mouse_vsd <- assay(mouse_vsd)
AP.mouse <- CoGAPS(mouse_vsd, nIterations=1000)
pNMF <-heatmap.2(as.matrix(AP.mouse),col=bluered, trace='none',
                 distfun=function(c) as.dist(1-cor(t(c))) ,
                 cexCol=1,cexRow=.5,scale = "row",
                 hclustfun=function(x) hclust(x, method="average"))

```