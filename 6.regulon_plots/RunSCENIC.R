#!/public2/Bio/pipeline/Toolkit_ycli/miniforge3/conda/envs/scenic/bin/Rscript
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

env <- tools::file_path_as_absolute(file.path(lib.dir, "scenic"))
.libPaths(file.path(env, "lib/R/library"))

args <- commandArgs(TRUE)

scenicOptions_file <- args[1]
exprMat_file <- args[2]
outdir <- args[3]

if (anyNA(c(scenicOptions_file, exprMat_file, outdir))) {
  warning("\n  Usage: ", env, "/bin/Rscript ", me, " <scenicOptions.Rds> <exprMat_filtered.Rds> <outdir> \n", call. = FALSE, immediate. = TRUE)
  quit()
}

library(SCENIC)
source(tools::file_path_as_absolute(file.path(lib.dir, "SCENIC_lib.R")), chdir = TRUE)

scenicOptions_file <- tools::file_path_as_absolute(scenicOptions_file)
exprMat_file <- tools::file_path_as_absolute(exprMat_file)

dir.create(outdir, FALSE, TRUE)
outdir <- tools::file_path_as_absolute(outdir)
setwd(outdir)

Message(">>>>> Load scenicOptions: ", scenicOptions_file)
scenicOptions <- readRDS(scenicOptions_file)

Message(">>>>> Load exprMat: ", exprMat_file)
exprMat <- readRDS(exprMat_file)

Message(">>>>> Prepare SCENIC")
Message("----> Check 'genesKept' again ")
genesKept <- loadInt(scenicOptions, "genesKept")
exprMat_filtered <- exprMat[genesKept, ]

Message("----> log2(exprMat_filtered + 1)")
exprMat_filtered_log <- as(log2(exprMat_filtered + 1), "CsparseMatrix")

Message(">>>>> Run SCENIC")
Message("----> Run SCENIC_1 coexNetwork2modules...")
runSCENIC_1_coexNetwork2modules(scenicOptions)
Message("----> Run SCENIC_2 createRegulons...")
minGenes <- scenicOptions@settings$createRegulons$minGenes %||% 20
runSCENIC_2_createRegulons(scenicOptions, minGenes = minGenes)
Message("----> Run SCENIC_3 scoreCells...")
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log, skipBinaryThresholds = FALSE, skipHeatmap = TRUE, skipTsne = TRUE)
Message("----> Run SCENIC_4 aucell_binarize...")
runSCENIC_4_aucell_binarize(scenicOptions, exprMat = exprMat_filtered_log, skipBoxplot = FALSE, skipHeatmaps = TRUE, skipTsne = TRUE)

Message(">>>>> Run tSNE")
tSNE_fileName <- runTSNE_AUC(scenicOptions, aucType = "AUC", onlyHighConf = FALSE) 
tSNE_bin_fileName <- runTSNE_AUC(scenicOptions, aucType = "Binary", onlyHighConf = FALSE)

Message(">>>>> Run UMAP")
UMAP_fileName <- runUMAP_AUC(scenicOptions, aucType = "AUC", onlyHighConf = FALSE)
UMAP_bin_fileName <- runUMAP_AUC(scenicOptions, aucType = "Binary", onlyHighConf = FALSE)

Message(">>>>> Tidy results")
motifEnrichment <- read.table(getOutName(scenicOptions,'s2_motifEnrichment'), sep="\t", header=T, stringsAsFactors=F)
motifEnrichment <- motifEnrichment[, c("motif", "motifDb", "nEnrGenes", "enrichedGenes", "highlightedTFs", "TFinDB", "TF_highConf", "TF_lowConf")]
colnames(motifEnrichment) <- c("Motif", "motifDB", "nEnrGenes", "enrichedGenes", "highlightedTFs", "TFinDB", "TF_highConf", "TF_lowConf")

SCENIC_results <- list(
  aucMat = AUCell::getAUC(loadInt(scenicOptions, 'aucell_regulonAUC')),
  binaryMat = loadInt(scenicOptions, 'aucell_binary_full'),
  cellInfo = loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo")),
  colVars = loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars")),
  RegulonGeneSet = loadInt(scenicOptions, 'aucell_regulons'),
  Regulon2Targets = getRegulon2Targets(scenicOptions),
  RegulonStat = statRegulon(scenicOptions),
  motifEnrichment = motifEnrichment,
  tf2target = getTF2Targets(scenicOptions),
  scenicOptions = scenicOptions,
  reductions = list(
    tSNE_AUC = readRDS(tSNE_fileName)$Y,
    tSNE_Binary = readRDS(tSNE_bin_fileName)$Y,
    UMAP_AUC = readRDS(UMAP_fileName),
    UMAP_Binary = readRDS(UMAP_bin_fileName)
  )
)

Message(">>>>> Save SCENIC results to: ", file.path(outdir, "SCENIC_final.Rds"))
saveRDS(SCENIC_results, file = file.path(outdir, "SCENIC_final.Rds"))

