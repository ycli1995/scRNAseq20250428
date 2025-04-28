
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

source("Seurat_lib.R", chdir=T)

library(future)
options(future.globals.maxSize = 100 * 1024 * 1024^2)
#plan("multiprocess", workers = 4)

obj = Load("obj.Rda")

load("markers.Rda")
param = yaml::yaml.load_file("RenameObject.yaml")

obj$seurat_clusters2 = factor(as.character(obj$seurat_clusters), levels = names(sort(table(obj$seurat_clusters), decreasing = T)))
obj@misc$color.cluster2 = obj@misc$color.cluster[levels(obj$seurat_clusters2)]

CalAvgExp(obj)

source("seuobj_lib.R", chdir = TRUE)
obj.markers$cluster = ReplaceEntries(obj.markers$cluster, param$rename$cluster$rename_map)
obj.markers$cluster = factor(as.character(obj.markers$cluster), levels = levels(obj$seurat_clusters2))

#obj.markers <- DoFindAllMarkers(obj)
#save(obj.markers, file = "markers.Rda")

if (FALSE) {
obj.markers$cluster = as.character(obj.markers$cluster)
#obj.markers$cluster = gsub(" |\\/", "_", obj.markers$cluster)
obj.markers$cluster = gsub("_", " ", obj.markers$cluster)
obj.markers$cluster = gsub(" S", "/S", obj.markers$cluster)
obj.markers$cluster = gsub(" E", "/E", obj.markers$cluster)
obj.markers$cluster = factor(obj.markers$cluster, levels(obj$seurat_clusters))
}
top <- FindTopMarker(obj.markers, top_num = 5, object = obj)

if (TRUE) {
PlotHeatmapPlot <- function(object, features = NULL, group.by = "seurat_clusters", is.use.name = TRUE, outfile = NULL, group.colors = NULL) {
		if ( is.null(group.colors) ) {
				group.colors <- switch(group.by,
									   "seurat_clusters" = object@misc$color.cluster,
									   "orig.ident" = object@misc$color.sample,
									   "Groups" = object@misc$color.group)
		} else {
				group.colors <- group.colors[levels(object@meta.data[[group.by]])]
		}
		p <- DoHeatmap( object = object, features = features, cells = NULL, group.by = group.by, group.colors = group.colors, combine = FALSE, raster = TRUE, label = FALSE)
		p <- p[[1]]
		p <- p + theme(legend.title = element_blank(), legend.text = element_text(size = 12, colour = "black"), axis.text.y.left = element_text(size = 12, colour = "black"))
			#p$layers[[2]] <- NULL
		if ( is.use.name ) {
				levels(p$data$Feature) <- FindFeaturesName(object, levels(p$data$Feature))
				levels(p$data$Feature) <- gsub("\\s(.*)", "", levels(p$data$Feature))
		}
		p <- p + scale_fill_gradient2(low = "#176BA1", mid = "#000000", high = "#FBAD3C")

		if ( is.null(outfile) ) {
			return(p)
		} else {
			h <- max(7, length(unique(features)) * 0.011 + 2.5 )
			w <- h * 4 / 3 + 0.5
			ggsave(p, file = outfile, width = w, height = h, limitsize = FALSE )
		}
}
}
features = unique(top$gene)
outpref = "Top"
group.by = "seurat_clusters2"
group.colors = obj@misc$color.cluster
is.use.name = TRUE
PlotHeatmapPlot(obj, features = features, group.by = group.by, outfile = paste0(outpref, ".Heatmap.pdf"), group.colors = group.colors, is.use.name = is.use.name)

