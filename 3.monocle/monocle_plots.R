library(monocle)
library(igraph)
library(dplyr)

add_lib = "Monocle_lib.R"
source(add_lib, chdir = T)

cls_col_name <- "seurat_clusters"
mnc_obj <- readRDS("Trajectory.obj.rds")
obj = Load("obj_used.Rda")
genes = unique(readLines("genes.list"))
BEAM_res = readRDS("BEAM.rds")

dir.create("add")
setwd("add")

library(Seurat)

name = "Trajectory."
show_branch_points = F
scale.size = 6
facet.scale <- 0.7

pData(mnc_obj)$Groups = obj@meta.data[colnames(mnc_obj), "Groups"]

## Trajectory #############################
mnc_obj2 = mnc_obj[, mnc_obj$Clusters %in% "EP1"]
mnc_obj2@reducedDimS = mnc_obj2@reducedDimS[, colnames(mnc_obj2)]
pData(mnc_obj2) = droplevels(pData(mnc_obj2))

p5 <- plot_cell_trajectory(mnc_obj2, color_by = "Clusters", show_branch_points = show_branch_points)  + scale_color_manual(values = obj@misc$color.cluster[levels(mnc_obj2$Clusters)])
ggsave( paste0( name, "Clusters.pdf" ), p5, width = scale.size, height = scale.size, limitsize = FALSE )

p5 <- plot_cell_trajectory(mnc_obj2, color_by = "Groups", show_branch_points = show_branch_points)  + scale_color_manual(values = obj@misc$color.group)
rc <- findRC( length(levels(pData(mnc_obj)$Groups)) )
p6 <- p5 + facet_wrap(~Groups, nrow = rc[1], scales = "free",)
ggsave( paste0( name, "Groups.pdf" ), p5, width = scale.size, height = scale.size, limitsize = FALSE );
ggsave( paste0( name, "Groups.facet.pdf" ),   p6, width = scale.size * rc[2] * facet.scale, height = scale.size * rc[1] * facet.scale, limitsize = FALSE )

p5 <- plot_cell_trajectory(mnc_obj2, color_by = "Samples", show_branch_points = show_branch_points)  + 
  scale_color_manual(values = obj@misc$color.sample) +
    guides(color = guide_legend(nrow=3,byrow=TRUE))
rc <- findRC( length(levels(pData(mnc_obj)$Samples)) )
p6 <- p5 + facet_wrap(~Samples, nrow = rc[1], scales = "free",)
ggsave( paste0( name, "Samples.pdf" ), p5, width = scale.size, height = scale.size, limitsize = FALSE );
ggsave( paste0( name, "Samples.facet.pdf" ),   p6, width = scale.size * rc[2] * facet.scale, height = scale.size * rc[1] * facet.scale, limitsize = FALSE )

## Gene expression
branch_point <- mnc_obj@auxOrderingData[[mnc_obj@dim_reduce_type]]$branch_points
branch_point_index <- FindBranchIndex(mnc_obj)
for (i in seq(branch_point)) {
  branch_name <- FindBranchName(ChangeBranchPointState(mnc_obj), branch_point[i])
  p4 <- plot_genes_branched_pseudotime(mnc_obj[genes, ], branch_point = i, color_by = "Groups", ncol = 4, branch_labels = branch_name) +
    scale_color_manual(values = obj@misc$color.group)
  ggsave(p4, file = paste0("Branch.", branch_point_index[branch_point[i]], ".genes_pseudotime.pdf"), width = 11, height = 9, limitsize = FALSE)
  for (g in genes) {
    p4 <- plot_genes_branched_pseudotime(mnc_obj[g, ], branch_point = i, color_by = "Groups", ncol = 1, branch_labels = branch_name) + 
      scale_color_manual(values = obj@misc$color.group)
    ggsave(p4, file = paste0("Branch.", branch_point_index[branch_point[i]], ".genes_pseudotime.", FindFeaturesName(obj, g), ".pdf"), width = 5, height = 4, limitsize = FALSE)
  }  
}



