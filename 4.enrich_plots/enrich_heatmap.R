
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)
library(circlize)

enrich_dir = "enrich"

all_cl = list.files(enrich_dir, pattern = "\\.glist$")
all_cl = gsub("\\.glist$", "", all_cl)

ref_path = read.table("ref.path.xls", header = T, sep = "\t", stringsAsFactors = F)
rownames(ref_path) = ref_path$Pathway

# KEGG
kegg_df = list()
for (i in all_cl) {
	kegg_df[[i]] = read.table(file.path(enrich_dir, "KO", paste0(i, ".bar_Gradient.xls")), header = T, sep = "\t", stringsAsFactors = F, check.names = F)
	#colnames(kegg_df[[i]])[4:5] = c("Cluster", "All")
	#kegg_df[[i]]$All = NULL
	kegg_df[[i]]$Cluster = i
}
kegg_df = do.call(rbind, kegg_df)
kegg_df$log10q = -log10(kegg_df$Qvalue)
top_df = kegg_df %>% group_by(Cluster) %>% top_n(5, log10q) %>% as.data.frame()

mat = dcast(kegg_df[kegg_df$Descrption %in% top_df$Descrption, ], Descrption ~ Cluster, value.var = "log10q")
rownames(mat) = mat[, 1]
mat[, 1] = NULL
mat = as.matrix(mat)[unique(top_df$Descrption), stringr::str_sort(all_cl, numeric = T)]
mat[is.na(mat)] = 0

col_fun = colorRamp2(c(0, min(floor(max(mat)), 20)), c("white", "red"))
hmp = Heatmap(
  mat,
  name = "-log10(Qvalue)",
  col = col_fun,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  show_column_names = T,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 8),
  column_title = "",
  column_title_side = "bottom",
  row_title = "",
  row_title_side = "left"
)
w = max(6, max(nchar(rownames(mat)) * 0.01 + 3.5))
h = 1.5 + nrow(mat) * 0.1
pdf("KEGG.Heatmap.pdf", width = w, height = h)
draw(hmp, merge_legend = T, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()

mat = dcast(kegg_df, Descrption ~ Cluster, value.var = "Pvalue")
rownames(mat) = mat[, 1]
mat[, 1] = NULL
mat = as.matrix(mat)[, stringr::str_sort(all_cl, numeric = T)]
#mat[is.na(mat)] = 1
mat2 = cbind(ref_path[rownames(mat), c(1:3)], as.data.frame(mat))
write.table(mat2, "KEGG_use.xls", sep = "\t", col.names = T, row.names = F, quote = F)

# GO
for (j in c("Cellular Component", "Biological Process", "Molecular Function")) {
  kegg_df = list()
  for (i in all_cl) {
    kegg_df[[i]] = data.table::fread(file.path(enrich_dir, "GO", paste0(i, ".bar_Gradient.xls")), header = T, sep = "\t", stringsAsFactors = F, check.names = F, data.table = F)
    kegg_df[[i]] = kegg_df[[i]][kegg_df[[i]]$class == j, ]
    kegg_df[[i]]$Cluster = i
  }
  kegg_df = do.call(rbind, kegg_df)
  kegg_df$log10q = -log10(kegg_df$Qvalue)
  top_df = kegg_df %>% group_by(Cluster) %>% top_n(5, log10q) %>% as.data.frame()
  mat = dcast(kegg_df[kegg_df$Descrption %in% top_df$Descrption, ], Descrption ~ Cluster, value.var = "log10q")
  rownames(mat) = mat[, 1]
  mat[, 1] = NULL
  mat = as.matrix(mat)[unique(top_df$Descrption), stringr::str_sort(all_cl, numeric = T)]
  mat[is.na(mat)] = 0

  col_fun = colorRamp2(c(0, min(floor(max(mat)), 50)), c("white", "red"))
  hmp = Heatmap(
    mat,        
    name = "-log10(Qvalue)",
    col = col_fun,
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = T,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    show_column_names = T,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 8),
    column_title = "",
    column_title_side = "bottom",
    row_title = "",
    row_title_side = "left"
  )
  w = max(8, max(nchar(rownames(mat)) * 0.01 + 4))
  h = 1.5 + nrow(mat) * 0.1
  j2 = gsub(" ", "_", j)
  outfile = paste0(j2, ".Heatmap.pdf")
  pdf(outfile, width = w, height = h)
  draw(hmp, merge_legend = T, padding = unit(c(2, 75, 2, 2), "mm"))
  dev.off()
  
  mat = dcast(kegg_df, Descrption ~ Cluster, value.var = "Pvalue")
  rownames(mat) = mat[, 1]
  mat[, 1] = NULL
  mat = as.matrix(mat)[, stringr::str_sort(all_cl, numeric = T)]
#  mat[is.na(mat)] = 0
  mat2 = cbind(data.frame(Description = rownames(mat), stringsAsFactors = F), as.data.frame(mat))
  write.table(mat2, paste0(j2, "_use.xls"), sep = "\t", col.names = T, row.names = F, quote = F)
}

