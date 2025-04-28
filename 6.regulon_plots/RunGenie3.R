#!/public2/Bio/pipeline/Toolkit_ycli/miniforge3/conda/envs/scenic/bin/Rscript
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

env <- tools::file_path_as_absolute(file.path(lib.dir, "scenic"))
.libPaths(file.path(env, "lib/R/library"))

args <- commandArgs(TRUE)

scenicOptions_file <- args[1]
exprMat_file <- args[2]
glist <- args[3]
index <- args[4]

if (anyNA(c(scenicOptions_file, exprMat_file, glist, index))) {
  warning("\n  Usage: ", env, "/bin/Rscript ", me, " <scenicOptions.Rds> <exprMat_filted.Rds> <glist> <index> \n", call. = FALSE, immediate. = TRUE)
  quit()
}

source(tools::file_path_as_absolute(file.path(lib.dir, "SCENIC_lib.R")), chdir = TRUE)

library(SCENIC)

scenicOptions_file <- tools::file_path_as_absolute(scenicOptions_file)
exprMat_file <- tools::file_path_as_absolute(exprMat_file)
glist <- tools::file_path_as_absolute(glist)

glist_dir <- dirname(glist)
glist_base <- basename(glist)
complete_file <- file.path(glist_dir, paste0("_", glist_base))

if (file.exists(complete_file)) {
  Message(">>>>> The run already complete: ", complete_file)
  quit()
}

Message(">>>>> Load scenicOptions: ", scenicOptions_file)
scenicOptions <- readRDS(scenicOptions_file)

Message(">>>>> Load exprMat: ", exprMat_file)
exprMat <- readRDS(exprMat_file)

Message(">>>>> Read gene list: ", glist)
genes.use <- readLines(con = glist)

Message(">>>>> RunGenie3")
exprMat_log <- as.matrix(log2(exprMat + 1))
RunGenie3(exprMat_log, scenicOptions, genes.use, index, nCores = 2)

write("complete", complete_file)

Message(">>>>> Done")

