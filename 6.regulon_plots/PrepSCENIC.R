#!/public2/Bio/pipeline/Toolkit_ycli/miniforge3/conda/envs/scenic/bin/Rscript
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

env <- tools::file_path_as_absolute(file.path(lib.dir, "scenic"))
.libPaths(file.path(env, "lib/R/library"))

args <- commandArgs(TRUE)

file <- args[1]
obj_file <- args[2]
outdir <- args[3]

if (anyNA(c(file, obj_file, outdir))) {
  warning("\n  Usage: ", env, "/bin/Rscript ", me, " <parameter.yaml> <obj_file> <outdir>\n", call. = FALSE, immediate. = TRUE) # nolint
  quit()
}

library(SCENIC)
source(tools::file_path_as_absolute(file.path(lib.dir, "SCENIC_lib.R")), chdir = TRUE)

file <- tools::file_path_as_absolute(file)
obj_file <- tools::file_path_as_absolute(obj_file)

handlers <- YAML_HANDLERS 
parameter <- yaml::yaml.load_file(file, handlers = handlers)

parameter$assay <- parameter$assay %||% "RNA"

parameter$scenicOptions$region <- parameter$scenicOptions$region %||% "both"
parameter$scenicOptions$genes.filter$minSamples.pct <- parameter$scenicOptions$genes.filter$minSamples.pct %||% 0.01
parameter$scenicOptions$genes.filter$minCountsPerGene.pct <- parameter$scenicOptions$genes.filter$minCountsPerGene.pct %||% 0.03

parameter$scenicOptions$RunGene3.threads <- parameter$scenicOptions$RunGene3.threads %||% 10
parameter$scenicOptions$createRegulons$minGenes <- parameter$scenicOptions$createRegulons$minGenes %||% 20
parameter$scenicOptions$main.threads <- parameter$scenicOptions$main.threads %||% 5

parameter$scenicOptions$group.by <- parameter$scenicOptions$group.by %||% c("Samples", "Groups", "Cluster")

print(parameter)

dir.create(outdir, FALSE, TRUE)
outdir <- tools::file_path_as_absolute(outdir)
setwd(outdir)

Message('>>>>> Loading Seurat object: ', obj_file)
obj <- Load(obj_file)
obj

Message(">>>>> 1. Initiate scenicOptions")
Message("----> Set options")
scenicOptions <- InitScenicOptions(parameter$scenicOptions, outdir = outdir)
scenicOptions@settings$createRegulons$minGenes <- parameter$scenicOptions$createRegulons$minGenes

Message("----> Set cellInfo")
cellInfo <- FetchData(obj, parameter$scenicOptions$group.by)
saveRDS(cellInfo, file = getDatasetInfo(scenicOptions, "cellInfo"))

Message("----> Set colVars")
colVars <- list()
for (i in colnames(cellInfo)) {
  colVars[[i]] <- CheckColorMap(cellInfo[, i], obj@misc$colors[[i]])
}
saveRDS(colVars, file = getDatasetInfo(scenicOptions, "colVars"))

Message(">>>>> 2. Fetch exprMat")
DefaultAssay(obj) <- parameter$assay
exprMat <- GetAssayData(obj, slot = "counts")

Message("----> Remove duplicated gene symbols")
gene_names <- FindFeaturesName(obj, rownames(exprMat))
gene_names <- setNames(gsub("\\s.*\\(.*", "", gene_names), names(gene_names)) ## 去除同名基因后缀的ID
gene_names <- gene_names[!duplicated(gene_names)]
exprMat <- exprMat[intersect(rownames(exprMat), names(gene_names)), , drop = FALSE]
gene_names <- gene_names[rownames(exprMat)]
rownames(exprMat) <- gene_names

saveRDS(exprMat, file = "exprMat.Rds")

Message(">>>>> 3. Filter gene")
parameter$scenicOptions$genes.filter$minSamples <- parameter$scenicOptions$genes.filter$minSamples %||% 
  parameter$scenicOptions$genes.filter$minSamples.pct * ncol(exprMat)
parameter$scenicOptions$genes.filter$minCountsPerGene <- parameter$scenicOptions$genes.filter$minCountsPerGene %||% 
  parameter$scenicOptions$genes.filter$minCountsPerGene.pct * ncol(exprMat)
message("minSamples: ", parameter$scenicOptions$genes.filter$minSamples)
message("minCountsPerGene: ", parameter$scenicOptions$genes.filter$minCountsPerGene)
genesKept <- geneFiltering(
  exprMat, 
  scenicOptions = scenicOptions,
  minCountsPerGene = parameter$scenicOptions$genes.filter$minCountsPerGene,
  minSamples = parameter$scenicOptions$genes.filter$minSamples
)
message(
  "Keep ", length(genesKept), " genes: \n",
  "  ", paste(head(genesKept), collapse = ", "), "..."
)
exprMat_filtered <- exprMat[genesKept, , drop = FALSE]
saveRDS(exprMat_filtered, file = "exprMat_filtered.Rds")

Message(">>>>> 4. Run Correlation")
exprMat_filtered <- as.matrix(exprMat_filtered)
gc(verbose = FALSE)
runCorrelation(exprMat_filtered, scenicOptions)

Message(">>>>> 5. Split genes")
genesSplit <- suppressWarnings(split(sort(rownames(exprMat_filtered)), 1:parameter$scenicOptions$RunGene3.threads))
for (i in 1:length(genesSplit)) {
  genes <- genesSplit[[i]]
  genes_file <- paste0("genesplit-", i)
  write.table(genes, file = genes_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

saveRDS(scenicOptions, file = "scenicOptions_raw.Rds")

Message(">>>>> Done")

