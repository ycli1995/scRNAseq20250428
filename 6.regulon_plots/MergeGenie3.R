#!/public2/Bio/pipeline/Toolkit_ycli/miniforge3/conda/envs/scenic/bin/Rscript
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

env <- tools::file_path_as_absolute(file.path(lib.dir, "scenic"))
.libPaths(file.path(env, "lib/R/library"))

args <- commandArgs(TRUE)

scenicOptions_file <- args[1]

if (anyNA(c(scenicOptions_file))) {
  warning("\n  Usage: ", env, "/bin/Rscript ", me, " <scenicOptions.Rds> \n", call. = FALSE, immediate. = TRUE)
  quit()
}

source(tools::file_path_as_absolute(file.path(lib.dir, "SCENIC_lib.R")), chdir = TRUE)

library(SCENIC)

scenicOptions_file <- tools::file_path_as_absolute(scenicOptions_file)

Message(">>>>> Load scenicOptions: ", scenicOptions_file)
scenicOptions <- readRDS(scenicOptions_file)

Message(">>>>> Merge Genie3")
MergeGenie3(scenicOptions)

Message(">>>>> Done")

