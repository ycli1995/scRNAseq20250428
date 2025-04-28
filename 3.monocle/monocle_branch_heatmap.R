
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

thres = 1e-7
cores = 1
use.q = TRUE

branch_point <- mnc_obj@auxOrderingData[[mnc_obj@dim_reduce_type]]$branch_points
branch_point_index <- FindBranchIndex(mnc_obj)

mats = list()
for (i in seq(branch_point)) {
  branch_name <- FindBranchName(ChangeBranchPointState(mnc_obj), branch_point[i])
  sig_gene_names <- select_sig_gene(BEAM_res[[i]], thres = thres, use.q = use.q)
  num_clusters <- min(length(sig_gene_names), 5)
  p5 <- plot_genes_branched_heatmap(mnc_obj[sig_gene_names, ], branch_point = i, num_clusters = num_clusters, cores = cores, show_rownames = T, return_heatmap = T, branch_labels = branch_name)
  mats[[i]] = list()
  mats[[i]]$mat = p5$heatmap_matrix
  mats[[i]]$row_anno = p5$annotation_row
  mats[[i]]$col_anno = p5$annotation_col
  mats[[i]]$hmcols = p5$hmcols
  mats[[i]]$col_anno_colors = p5$annotation_colors
}
saveRDS(mats, file = "monocle_branch_heatmap_mats.Rds")

library(ComplexHeatmap)
genes = FindFeaturesName(obj, genes)
heatmap_legend_param <- list(title_position = 'leftcenter-rot', legend_height = unit(3, 'cm'))

for (i in seq(mats)) {
  col_anno = mats[[i]]$col_anno
  col_anno$`Cell Type` = factor(col_anno$`Cell Type`, levels = unique(col_anno$`Cell Type`))
  col_anno_colors = mats[[i]]$col_anno_colors
  col_anno_colors$`Cell Type` = col_anno_colors$`Cell Type`[unique(col_anno$`Cell Type`)]
  top_anno = HeatmapAnnotation(df = col_anno, col = col_anno_colors)
  
  left_colors = scales::hue_pal()(length(levels(mats[[i]]$row_anno$Cluster))) %>%
    setNames(levels(mats[[i]]$row_anno$Cluster))
  left_anno = HeatmapAnnotation(df = mats[[i]]$row_anno, col = list(Cluster = left_colors), which = "row")
  ## Get label genes
  g_used = intersect(rownames(mats[[i]]$mat), genes)
  anno = anno_mark(at = match(g_used, rownames(mats[[i]]$mat)), labels = g_used, which = "row")
  
  col_fun = circlize::colorRamp2(seq(min(mats[[i]]$mat), max(mats[[i]]$mat), length.out = length(mats[[i]]$hmcols)), mats[[i]]$hmcols)
  hmp = Heatmap(
    mats[[i]]$mat, col = col_fun, name = "Expression",
    top_annotation = top_anno,
    left_annotation = left_anno,
    right_annotation = rowAnnotation(mark = anno),
    row_split = mats[[i]]$row_anno$Cluster,
    column_split = c(rep(1, 100), rep(2, 100)),
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    row_title = NULL,
    column_title = NULL
    #heatmap_legend_param = heatmap_legend_param   
  )

  pdf(paste0("Branch.", i, ".Heatmap.pdf"), height = 6, width = 8)
  draw(hmp)
  dev.off()
  svglite::svglite(paste0("Branch.", i, ".Heatmap.svg"), height = 6, width = 8)
  draw(hmp)
  dev.off()
}






