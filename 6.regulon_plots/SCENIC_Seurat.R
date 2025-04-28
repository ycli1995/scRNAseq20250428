#!/public2/Bio/pipeline/Toolkit_ycli/miniforge3/conda/envs/scenic/bin/Rscript
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

env <- tools::file_path_as_absolute(file.path(lib.dir, "scenic"))
.libPaths(file.path(env, "lib/R/library"))

args <- commandArgs(TRUE)

obj_file <- args[1]
scenic_file <- args[2]
outdir <- args[3]

if (anyNA(c(obj_file, scenic_file, outdir))) {
  warning("\n  Usage: ", env, "/bin/Rscript ", me, " <obj_file> <SCENIC_final.Rds> <outdir>\n", call. = FALSE, immediate. = TRUE)
  quit()
}

library(SCENIC)
#source(tools::file_path_as_absolute(file.path(lib.dir, "Heatmap_lib.R")), chdir = TRUE)
source(tools::file_path_as_absolute(file.path(lib.dir, "SCENIC_lib.R")), chdir = TRUE)

obj_file <- tools::file_path_as_absolute(obj_file)
scenic_file <- tools::file_path_as_absolute(scenic_file)

dir.create(outdir, FALSE, TRUE)
outdir <- tools::file_path_as_absolute(outdir)
setwd(outdir)

Message('>>>>> Loading Seurat object: ', obj_file)
obj <- Load(obj_file)
obj

Message('>>>>> Loading SCENIC results: ', scenic_file)
SCENIC_results <- readRDS(scenic_file)

Message('>>>>> 1. SCENIC tables')
WriteTable(cbind(id = rownames(SCENIC_results$aucMat), SCENIC_results$aucMat), "RegulonActivity.xls")
WriteTable(cbind(id = rownames(SCENIC_results$binaryMat), SCENIC_results$binaryMat), "BinaryRegulonActivity.xls")

aucMat_avg <- list()
binaryMat_avg <- list()

group_by <- colnames(SCENIC_results$cellInfo)
# For report
write(group_by, "groupby.list")
write(levels(as.factor(obj$Samples)), "sample.list")

for (i in group_by) {
  aucMat_avg[[i]] <- CalAvgExp(SCENIC_results$aucMat, group.by = as.factor(SCENIC_results$cellInfo[, i]), expm1 = FALSE)
  WriteTable(cbind(id = rownames(aucMat_avg[[i]]), aucMat_avg[[i]]), paste0("RegulonActivity.avg.", i, ".xls"))
  binaryMat_avg[[i]] <- CalAvgExp(SCENIC_results$binaryMat, group.by = as.factor(SCENIC_results$cellInfo[, i]), expm1 = FALSE)
  WriteTable(cbind(id = rownames(binaryMat_avg[[i]]), binaryMat_avg[[i]]), paste0("BinaryRegulonActivity.avg.", i, ".xls"))
}

WriteTable(SCENIC_results$RegulonStat, "RegulonStat.xls")
WriteTable(SCENIC_results$motifEnrichment, "motifEnrichment.xls")
WriteTable(SCENIC_results$tf2target, "tf2target.xls")

RegulonGeneSet <- data.frame(
  regulon = names(SCENIC_results$RegulonGeneSet),
  target = unlist(lapply(SCENIC_results$RegulonGeneSet, paste, collapse = ","))							 
)
WriteTable(RegulonGeneSet, "RegulonGeneSet.xls")

Message('>>>>> 2. SCENIC reductions')
## For online
dim_df <- FetchData2(obj, group_by)

for (i in names(SCENIC_results$reductions)) {
  Message('----> ', i)
  colnames(SCENIC_results$reductions[[i]]) <- paste0(i, "_", 1:ncol(SCENIC_results$reductions[[i]]))
  dim_df <- cbind(dim_df, SCENIC_results$reductions[[i]][Cells(obj), ])
  obj[[i]] <- CreateDimReducObject(SCENIC_results$reductions[[i]][Cells(obj), ], assay = DefaultAssay(obj))
  for (j in group_by) {
    params <- list()
    obj@meta.data[, j] <- SCENIC_results$cellInfo[Cells(obj), j]
    params[[j]][["cols"]] <- SCENIC_results$cellVars[[j]]
    PlotCluster(obj, reduction = i, group.by = j, params = params, outpref = paste0(i, ".", j), combine = TRUE, raster = FALSE)
  }
}
for (i in c("tsne", "umap")) {
  if (i %in% Reductions(obj)) {
    dim_df <- cbind(dim_df, obj[[i]]@cell.embeddings)
  }
}
WriteTable(dim_df, "SCENIC_reductions.tmp")

Message('>>>>> 3. Add TF and target annotations')
fdata <- obj@misc$fdata
fdata$dash <- gsub("_", "-", rownames(fdata))
fdata$Regulon <- ""
if ("gene_symbols" %in% colnames(fdata)) {
  fdata$name = fdata$gene_symbols
}
for (i in names(SCENIC_results$RegulonGeneSet)) {
  fdata$Regulon[fdata$name %in% SCENIC_results$RegulonGeneSet[[i]]] <- paste0(fdata$Regulon[fdata$name %in% SCENIC_results$RegulonGeneSet[[i]]], i, "; ")
}
fdata$Regulon[fdata$Regulon %in% ""] <- "--"
fdata$Type <- ""
fdata$Type[fdata$name %in% SCENIC_results$tf2target$tf] <- "TF;"
fdata$Type[fdata$name %in% SCENIC_results$tf2target$target] <- paste0(fdata$Type[fdata$name %in% SCENIC_results$tf2target$target], "Target;")
fdata$Type[fdata$Type %in% ""] <- "--"

fdata <- cbind(GeneID = rownames(fdata), GeneName = fdata$name, fdata[, c("Type", "Regulon", "dash")])
rownames(fdata) <- fdata$dash
fdata$dash <- NULL

for (i in group_by) {
  mean_exp <- CalAvgExp(obj, group.by = i)
  mean_exp <- cbind(fdata[rownames(mean_exp), ], mean_exp)
  WriteTable(mean_exp, paste0("AllGene.avg.", i, ".xls"))
  mat <- GetAssayData(obj)
  if (inherits(mat, "dgCMatrix")) {
    mat@x[mat@x > 0] <- 1
    mat <- drop0(mat)
  } else {
    mat[mat > 0] <- 1
  }
  #mat[mat > 0] <- 1
  g <- SCENIC_results$cellInfo[, i]
  pct_exp <- CalAvgExp(mat, group.by = g, expm = FALSE)
  pct_exp <- cbind(fdata[rownames(pct_exp), ], pct_exp)
  WriteTable(pct_exp, paste0("AllGene.pct.", i, ".xls"))
}

Message('>>>>> 4. Tables for network')
node_annot <- mean_exp[, 1:3]
node_annot <- node_annot[!node_annot$Type %in% "--", , drop = FALSE]
WriteTable(node_annot, "Network.Node.Annot.xls")

### Old version pipeline
tf2target <- SCENIC_results$tf2target
edge <- tf2target[, 2:4]
WriteTable(edge, "tf_target.egde.tsv")

node <- data.frame(node = unique(c(edge$tf, edge$target)), stringsAsFactors = FALSE)
node$annot <- "target"
node$annot[node$node %in% edge$tf] <- "tf"
WriteTable(node, "tf_target.node.tsv")

Message('>>>>> 5. Get SCENIC Seurat obj')
tf_genes <- gsub("_extend.*", "", rownames(SCENIC_results$aucMat))
tf_genes <- gsub(" .*", "", tf_genes)
tf_genes <- FindFeaturesID(obj, tf_genes)
#obj <- DietSeurat(obj, dimreducs = Reductions(obj))
for (i in SeuratObject::Assays(obj)) {
  if (all(tf_genes %in% rownames(obj[[i]]))) {
    obj[[i]] <- subset(obj[[i]], features = tf_genes)
  }
}
obj[['SCENIC']] <- CreateAssayObject(data = SCENIC_results$aucMat)
rownames(obj[['SCENIC']]@data) <- gsub("\\-extended", "_extended", rownames(obj[['SCENIC']]@data))
obj[['SCENIC']]@counts <- SCENIC_results$binaryMat
obj <- DietSeurat(obj, assays = Assays(obj), dimreducs = Reductions(obj))
obj@misc$counts <- NULL
save(obj, file = "SCENIC_obj.Rda")

Message('>>>>> Done')

