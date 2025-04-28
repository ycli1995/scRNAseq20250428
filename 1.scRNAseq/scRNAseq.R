
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

source("Seurat_lib.R", chdir = T)

.PlotClusterStat = function(object, stat.what = "seurat_clusters", group.by = "orig.ident", color.st = NULL, color.gb = NULL, outpref = NULL, ...){
		if ( class(object) == "Seurat" ) {
				metadata <- object@meta.data
		} else {
				metadata <- object
		}
		if ( is.null(color.st) ) {
				if ( "misc" %in% slotNames(object) && exists(stat.what, object@misc) ) {
						color.st <- object@misc[[stat.what]]
				} else {
				color.st <- switch(stat.what,
						"seurat_clusters" = object@misc$color.cluster,
						"orig.ident" = object@misc$color.sample,
						"Groups" = object@misc$color.group
				)
				}
		}
		if ( is.null(color.gb) ) {
				if ( "misc" %in% slotNames(object) && exists(group.by, object@misc) ) {
						color.gb <- object@misc[[group.by]]
				} else {
				color.gb <- switch(group.by, 
						"seurat_clusters" = object@misc$color.cluster,
						"orig.ident" = object@misc$color.sample,
						"Groups" = object@misc$color.group
				)
				}
		}
		name.st <- switch(stat.what, "seurat_clusters" = "Cluster", "orig.ident" = "Samples", stat.what)
		name.gb <- switch(group.by,  "seurat_clusters" = "Cluster", "orig.ident" = "Samples", group.by)
		stat.what <- as.name(stat.what)
		group.by  <- as.name(group.by)

		stat_sample <- metadata %>%
				group_by(!! name.gb := !! group.by, !! name.st := !! stat.what) %>%
				summarise("Number of cells" = n())

		p <- list()
		p[["by"]] <- ggplot(stat_sample, aes_(x = as.name(name.gb), y = ~ `Number of cells`, fill = as.name(name.st)))
		p[["in"]] <- ggplot(stat_sample, aes_(x = as.name(name.st), y = ~ `Number of cells`, fill = as.name(name.gb)))
		if ( ! is.null(color.st) ) p[["by"]] <- p[["by"]] + scale_fill_manual(values = color.st)
		if ( ! is.null(color.gb) ) p[["in"]] <- p[["in"]] + scale_fill_manual(values = color.gb)
		geom_stack <- geom_bar(stat = "identity", position = 'stack')
		geom_fill  <- geom_bar(stat = "identity", position = "fill" )

		if ( is.null(outpref) ) {
				outpref <- paste0( name.st, ".stat")
		}
		for ( i in names(p) ) {
				p[[i]] <- p[[i]] + bar_theme_default()
				if (i == "by") w = length(unique(stat_sample[, name.gb])) * 1.75 + 2.5
				if (i == "in") w = length(unique(stat_sample[, name.st])) * 1.75 + 2.5
				ggsave( p[[i]] + geom_stack, file = paste0( outpref, ".", i, name.gb, ".pdf"), height = 6, width = w )
				ggsave( p[[i]] + geom_fill + ylab("Fraction of Cells"),  file = paste0( outpref, ".", i, name.gb, ".pct.pdf"), height = 6, width = w )
		}
}

StatCluster = function(object, group.by = "orig.ident", outpref = "Cluster.stat", stat.what = "seurat_clusters", assay = DefaultAssay(object), ...){
		.StatCluster(object, stat.what = stat.what, outpref = outpref, assay = assay)
		.StatCluster_by(object, stat.what = stat.what, group.by = group.by, outpref = outpref)
		.PlotClusterStat(object, stat.what = stat.what, group.by = group.by, outpref = outpref, ...)
}

obj = Load("obj_renamed.Rda")
genes = read.table("genes.list", sep = "\t", stringsAsFactors = F)$V1

# 1. umap
obj$seurat_clusters2 = factor(as.character(obj$seurat_clusters), levels = names(sort(table(obj$seurat_clusters), decreasing = T)))
obj@misc$color.cluster2 = obj@misc$color.cluster[levels(obj$seurat_clusters2)]

PlotCluster(obj, outpref = "UMAP", p2.group.by = "seurat_clusters2", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, p2.label = F, outpref = "UMAP.nolabel", p2.group.by = "seurat_clusters2", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, outpref = "UMAP.Group", p1.group.by = "Groups", p2.group.by = "seurat_clusters2", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, p2.label = F, outpref = "UMAP.nolabel.Group", p1.group.by = "Groups", p2.group.by = "seurat_clusters2", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, "tsne", outpref = "TSNE", p2.group.by = "seurat_clusters2", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, "tsne", p2.label = F, outpref = "TSNE.nolabel", p2.group.by = "seurat_clusters2", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, "tsne", outpref = "TSNE.Group", p1.group.by = "Groups", p2.group.by = "seurat_clusters2", p2.color = obj@misc$color.cluster2)
PlotCluster(obj, "tsne", p2.label = F, outpref = "TSNE.nolabel.Group", p1.group.by = "Groups", p2.group.by = "seurat_clusters2", p2.color = obj@misc$color.cluster2)

# 2. Stat
StatCluster(obj, stat.what = "seurat_clusters2", color.st = obj@misc$color.cluster2)
StatCluster(obj, "Groups", stat.what = "seurat_clusters2", color.st = obj@misc$color.cluster2)

#obj2 = subset(obj, orig.ident %in% c("Pre-Cart-T-2", "Pre-Cart-T-3"))
#obj2@meta.data = droplevels(obj2@meta.data)
#obj2@misc$color.sample = obj2@misc$color.sample[levels(obj2$orig.ident)]

#StatCluster(obj2, outpref = "Cluster.subset.stat")

# 3. DotPlot
p = PlotDotPlot(obj, FindFeaturesID(obj, genes)) + 
  scale_colour_gradientn(colors = RColorBrewer::brewer.pal(9, 'OrRd')) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 0)
  ) + 
  theme(legend.box = "horizontal")
p$data$id = factor(as.character(p$data$id), levels = rev(levels(p$data$id))) 
p$data$features.plot = factor(as.character(p$data$features.plot), levels = rev(levels(p$data$features.plot)))

w <- max(7, ceiling(length(genes)) * 0.35 + 2)
h = max(2.75, length(unique(obj@meta.data[["seurat_clusters"]])) * 0.4)
ggsave(p, file = "DotPlot.pdf", width = w, height = h, limitsize = FALSE )

WriteTable(p$data, "DotPlot.xls")
