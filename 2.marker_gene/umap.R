
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

source("Seurat_lib.R", chdir=T)

obj = Load("obj.Rda")

features = read.table("genes.list", header=F, stringsAsFactors=F)$V1
features = FindFeaturesID(obj, features)
#features = features[features != "FBgn0011676"]

PlotFeaturePlot(obj, features, reduction = "umap", is.combine = FALSE, cols = c("#B1DAE9", "#B6212D"))

#features = read.table("../9.umap/gene.list", header=F, stringsAsFactors=F)$V1

#setwd("../10.violin")

#table(features %in% obj@misc$fdata$merge_name)

#fdata = obj@misc$fdata[rownames(obj[['RNA']]@data), ]

#rownames(obj[['RNA']]@data) = fdata$merge_name
#p1 = VlnPlot(obj, features = features, group.by = "seurat_clusters", stack = TRUE, cols = obj@misc$color.cluster, fill.by = "ident")
#p1$data$ident = factor(p1$data$ident, levels = rev(levels(p1$data$ident)))
#p1 = p1 + NoLegend() + 
#  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 10, angle = 90, face = "plain", hjust = 0, vjust = 0.5), axis.title.y = element_blank())
#ggsave(p1, file = "Violin.pdf", width = 16, height = 4)


#obj$seurat_clusters = factor(obj$seurat_clusters, levels = rev(levels(obj$seurat_clusters)))
#PlotDotPlot(obj, rev(features), outfile = "DotPlot.pdf", cols = c("#B1DAE9", "darkred"))
