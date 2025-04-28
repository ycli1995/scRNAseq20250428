
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

source("Seurat_lib.R", chdir = T)

k = 10
obj = Load("obj_used.Rda")
#CNV_res = readRDS("0.infercnv.All.statCNV.rds")
gene_order = read.table("gene_order.txt", header = F, sep = "\t", stringsAsFactors = F, row.names = 1)
colnames(gene_order) = c("chr", "start", "stop")
gene_order$chr = factor(as.character(gene_order$chr), levels = unique(gene_order$chr))

## CNV score
CNV = readRDS("run.final.infercnv_obj")
mat = CNV@expr.data[, Cells(obj)]
mat = mat[intersect(rownames(gene_order), rownames(mat)), ]

rowData = gene_order[rownames(mat), "chr", drop = F]

obj$CNV_score = (mat - 1) ^ 2 %>% Matrix::colMeans()

obj$infercnv_grouping = "Obs" 
obj$infercnv_grouping[obj$orig.ident %in% "Normal-Ovary-1"] = "Ref" 

## Kmeans
obj$infercnv_kmeans = "Ref"

km = kmeans(t(mat[, Cells(obj)[obj$infercnv_grouping %in% "Obs"]]), k)$cluster

obj$infercnv_kmeans[obj$infercnv_grouping %in% "Obs"] = as.character(km)
obj$infercnv_kmeans = factor(as.character(obj$infercnv_kmeans), levels = c("Ref", as.character(1:k)))
obj@misc$infercnv_kmeans = SetColor(obj$infercnv_kmeans, "tsne", "set1")

PlotVlnPlot(obj, "CNV_score", outpref = "Violin", group.by = "infercnv_kmeans", cols = obj@misc$infercnv_kmeans)

StatCluster(obj, group.by = "orig.ident", outpref = "Kmeans.stat", stat.what = "infercnv_kmeans", color.st = obj@misc$infercnv_kmeans)

## Prep for heatmap
colData = obj@meta.data
colData$Samples = colData$orig.ident

colors = list(Samples = obj@misc$color.sample, Groups = obj@misc$color.group, infercnv_kmeans = obj@misc$infercnv_kmeans, chr = SetColor(rowData$chr, "onlinereport"))

all_results = list(CNAmat = mat, colData = colData, rowData = rowData, colors = colors)
saveRDS(all_results, file = "infercnv.Rds")




