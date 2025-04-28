## Welcome GD_inferCNV.R ##
## Version: v1.0 --- lhsong@genedenovo.com ##
inferCNV_version = 1.0
source('Start.r')

## Dell arguments ##
args = commandArgs()
me   = normalizePath(sub('--file=', '', grep('--file=', args, value = T)), '/')
args = args[-seq(grep('--args', args))]

file = args[1]
obj  = args[2]
name = args[3]
main = 'Usage v1.0: /Bio/Bin/pipeline/InferCNV/dev/Rsc4 GD_inferCNV.R <inferCNV_paramter.yml> <obj.Rda> <name>\n Introduction: InferCNV analysis.'
if (is.na(file)) stop(Er(main))
if (is.na(name)) name = 'infercnv'

message(Wa('==>', timer(), 'Run: ', me, '<=='))
message(Sa('-->', timer(), 'Welcome to use GD_inferCNV', Pa('version:', inferCNV_version), '! <--'))
reticulate::use_python('/Bio/Bin/pipeline/InferCNV/dev/python', T)
file = normalizePath(file, '/', T)
obj  = normalizePath(obj, '/', T)

# lib #
suppressMessages(library(infercnv)) ; suppressMessages(library(Seurat))
suppressMessages(library(R.utils))  ; suppressMessages(library(future.apply))
suppressMessages(library(Matrix))   ; suppressMessages(library(ggplot2))
suppressMessages(library(reshape2)) ; suppressMessages(library(patchwork))
suppressMessages(library(Signac))

## parameter ##
message(Sa('-->', timer(), 'read config:', Pa(file), '<--'))
parameter  = yaml.load_file(file, handlers = handlers)
gene_order = normalizePath(parameter$gene_order, '/', T)
col_use    = parameter$colname_use  %||% 'orig.ident'
plot_cols  = parameter$plot_cols    %||% c('orig.ident', 'seurat_clusters', 'subcluster', 'kmeans')
ref_group  = parameter$ref_group
mode       = parameter$mode         %||% 'subcluster'
method     = parameter$mode_method  %||% 'leiden'
pvalue     = parameter$mode_pvalue  %||% 1e-6
threads    = parameter$threads      %||% 10
kmeans     = parameter$kmeans       %||% 5
topK       = parameter$topK         %||% 1
cut_off    = parameter$cut_off      %||% .1
window     = parameter$window       %||% 101
outdir     = parameter$outdir       %||% dirname(file)
# check essential para #
mode       = as.character(mode[1])
method     = as.character(method[1])
col_use    = as.character(col_use[1])
methods    = c('leiden', 'random_trees', 'qnorm', 'pheight', 'qgamma', 'shc')
if (!sum(mode %in% c('sample', 'subcluster'))) 
  stop(Er('!!!', timer(), 'run mode must be', Pa('[ sample | subcluster ]'), 'but input mode is:', Pa(mode), '!!!'))
if (mode %in% 'subcluster' & !sum(method %in% methods)) 
  stop(Er('!!!', timer(), 'subcluster method must be', Pa('[', paste(methods, collapse = ' | '), ']'), 'but input method is:', Pa(method), '!!!'))
# function #
cores      = function() threads
plan(multicore(workers = cores()))

# dir #
message(Sa('-->', timer(), 'Process outdir (make sure you have permission):', Pa(outdir), '<--'))
dir.create(outdir, F, T); setwd(outdir)
system('rm -f GD_inferCNV.done ; touch GD_inferCNV.undone')
wdir = getwd()

# 0. Load obj #
message(Sa('-->', timer(), '0. Read Seurat object (.Rda):', Pa(obj), 'for sample:', Pa(name), '<--'))
load(obj)
# check ref_group
outer = setdiff(ref_group, obj@meta.data[[col_use]])
if (length(outer)) {
  message(No('-!!', timer(), 'End: there is no ref_group:', Pa(outer, collapse = ','), 'in obj@meta.data', Pa(col_use), '!!-'))
  system('rm -f GD_inferCNV.undone ; touch GD_inferCNV.done') ; q('no')
}
# initail subcluster and kmeans
if ('subcluster' %in% colnames(obj@meta.data))
  colnames(obj@meta.data)[colnames(obj@meta.data) %in% 'subcluster'] = 'subcluster-bac'
obj$subcluster = NA
if ('kmeans' %in% colnames(obj@meta.data))
  colnames(obj@meta.data)[colnames(obj@meta.data) %in% 'kmeans'] = 'kmeans-bac'
obj$kmeans = NA
# check plot_cols
outer     = setdiff(plot_cols, colnames(obj@meta.data))
plot_cols = intersect(plot_cols, colnames(obj@meta.data))
if (length(outer))
  message(No('--!', timer(), 'There is no plot_cols:', Pa(outer, collapse = ','), 'in obj@meta.data, ignored !--'))
if (!length(plot_cols)) 
  stop(Er('!!!', timer(), 'No availble plot_cols in obj@meta.data !!!'))

# 1. Run infercnv #
message(Sa('-->', timer(), '1. Run infercnv in mode:', Pa(mode), if(mode %in% 'subcluster') Pa(method), '<--'))
tfile = paste0(name, '.cellAnnotations.txt')
annot = paste0(name, '-', ifelse(obj@meta.data[[col_use]] %in% ref_group, 'Ref', 'Obs'))
write.table(data.frame(cell = Cells(obj), annot = annot), tfile, sep = '\t', row.names = F, col.names = F)
message(Sa('-->', timer(), 'saved', Pa(tfile), 'at', Pa(wdir), '<--'))
raw   = paste0('raw_', name)
dir.create(raw, F)
# main run
message(Sa('-->', timer(), 'infercnv analysing ... (may take a long time) <--'))
rfile = paste0('0.infercnv.', mode, if(mode %in% 'subcluster') paste0('.', method), '.', name, '.rds')
if (!file.exists(rfile)) { 
  counts = if('RNA' %in% names(obj@assays)) obj@assays$RNA@counts else obj@assays$Spatial@counts
  ifcObj = CreateInfercnvObject(counts, annotations_file = tfile, gene_order_file = gene_order, ref_group_names = paste0(name, '-Ref'), chr_exclude = NULL)
  rm(counts)
  ifcObj = infercnv::run(ifcObj, cutoff = cut_off, analysis_mode = mode, 
                         HMM = T, HMM_type = 'i6', window_length = window,
                         tumor_subcluster_partition_method = method, tumor_subcluster_pval = pvalue,
                         plot_probabilities = F, plot_steps = F, out_dir = raw,
                         num_threads = threads, output_format = 'pdf', cluster_references = F, cluster_by_groups = T, denoise = T)
  add_to_seurat(, raw)
  saveRDS(ifcObj, rfile)
  message(Sa('-->', timer(), 'saved', Pa(rfile), 'at', Pa(wdir), '<--'))
  # mv obs and ref
  rn = file.rename(paste0(raw, '/infercnv.', c('observations.txt'  , 'references.txt')),
                   paste0(raw, '/0.infercnv.', name, '.', c('observations.txt', 'references.txt')) )
  # gzip results
  gz = future_lapply(list.files(raw, '.txt$', recursive = T, full.names = T), function(f) 
    gzip(f, overwrite = T) )
} else {
  message(Sa('-->', timer(), 'read from previous results:', Pa(rfile), '<--'))
  ifcObj = readRDS(rfile)
}
# gene state
tfile = paste0('1.infercnv.', name, '.annot.xls')
if (!file.exists(tfile)) {
  HMM   = read.table(list.files(raw, 'HMM_CNV_predictions.*.pred_cnv_genes.*', full.names = T)[1], sep = '\t', header = T, quote = '"')
  write.table(HMM, tfile, sep = '\t', row.names = F, quote = F)
  message(Sa('-->', timer(), 'saved', Pa(tfile), 'at', Pa(wdir), '<--'))
} else
  message(Sa('-->', timer(), 'checked previous HMM annot:', Pa(tfile), '<--'))
# chr state
tfile = paste0('1.infercnv.', name, '.chr.xls')
if (!file.exists(tfile)) {
  Chr   = read.table(paste0(raw, '/map_metadata_from_infercnv.txt.gz'), sep = '\t', header = T, quote = '"')
  write.table(cbind(barcode = rownames(Chr), Chr), tfile, sep = '\t', row.names = F, quote = F)
  message(Sa('-->', timer(), 'saved', Pa(tfile), 'at', Pa(wdir), '<--'))
} else
  message(Sa('-->', timer(), 'checked previous chr annot:', Pa(tfile), '<--'))
# subcluster
subcluster = do.call(rbind, lapply(ifcObj@tumor_subclusters$subclusters, function(clust) 
  do.call(rbind, lapply(seq(clust), function(i) 
    data.frame(Barcode = names(clust[[i]]), SubClust = names(clust[i])) ))))
# make ref subcluster -> reference #
levels = unlist(lapply(ifcObj@tumor_subclusters$subclusters, names))
obj$subcluster[match(subcluster$Barcode, Cells(obj))]   = subcluster$SubClust
obj$subcluster[obj@meta.data[[col_use]] %in% ref_group] = 'Ref'
levels[levels %in% setdiff(levels, obj$subcluster)]     = 'Ref'
obj$subcluster = factor(obj$subcluster, unique(levels))
# kmeans for predict
message(Sa('-->', timer(), 'kmeans ... <--'))
kmean = factor(kmeans(t(ifcObj@expr.data), kmeans)$cluster)
pred  = data.frame(observation = as.numeric(table(kmean[unlist(ifcObj@observation_grouped_cell_indices)])),
                   reference   = as.numeric(table(kmean[unlist(ifcObj@reference_grouped_cell_indices  )])) )
pred  = cbind(kmean = 1:nrow(pred), pred[order(-pred$reference),])
p1    = ggplot(melt(pred, id.vars = 'kmean'), aes(kmean, value, fill = variable)) + 
        geom_col() + labs(x = 'Kmeans', y = 'Frequency', fill = '') +
        scale_y_continuous(expand = expansion(c(.01, .05))) + theme_classic() +
        theme(text = element_text(size = 13))
p2    = ggplot(melt(pred, id.vars = 'kmean'), aes(kmean, value, fill = variable)) +  
        geom_col(position = 'fill') + labs(x = 'Kmeans', y = 'Proportion', fill = '') +
        scale_y_continuous(expand = expansion(c(.01, .05))) + theme_classic() +
        theme(text = element_text(size = 13))
p     = p1 + p2
pfile = paste0('1.infercnv.', name, '.predict.pdf')
ggsave(pfile, p, width = 14, height = 6, limitsize = F)
message(Sa('-->', timer(), 'saved', Pa(pfile), 'at', Pa(getwd()), '<--'))
rm(p1, p2, p)
pred$Total       = pred$observation + pred$reference
pred$observation = paste0(pred$observation, ' (', makePct(pred$observation / pred$Total),    ')')
pred$reference   = paste0(pred$reference  , ' (', makePct(pred$reference   / pred$Total),    ')')
pred$Total       = paste0(pred$Total      , ' (', makePct(pred$Total       / sum(pred$Total)), ')')
tfile = paste0('1.infercnv.', name, '.kmeans.xls')
write.table(pred, tfile, sep = '\t', row.names = F, quote = F)
message(Sa('-->', timer(), 'saved', Pa(tfile), 'at', Pa(wdir), '<--'))
# order kmean
levels(kmean)[as.numeric(rownames(pred))] = pred$kmean
kmean = as.character(kmean)
kmean[unlist(ifcObj@reference_grouped_cell_indices)] = '1'
obj$kmeans[match(colnames(ifcObj@expr.data), Cells(obj))] = kmean
# cnv predict
CNV   = data.frame(Cells = colnames(ifcObj@expr.data), Kmeans = kmean, Subcluster = 'none', Prediction = 'aneuploid')
CNV$Prediction[ kmean %in% seq(topK) ] = 'diploid'
CNV$Subcluster = subcluster$SubClust[match(CNV$Cells, subcluster$Barcode)]
tfile = paste0('1.infercnv.', name, '.predict.xls')
write.table(CNV, tfile, sep = '\t', row.names = F, quote = F)
message(Sa('-->', timer(), 'saved', Pa(tfile), 'at', Pa(wdir), '<--'))
rm(kmean, pred)
# stat cell type
stat           = table(CNV$Subcluster, CNV$Prediction)
stat           = data.frame(rbind(Total = apply(stat, 2, sum), stat))
stat$Total     = stat$aneuploid + stat$diploid
stat$aneuploid = paste0(stat$aneuploid, ' (', makePct(stat$aneuploid / stat$Total),    ')')
stat$diploid   = paste0(stat$diploid  , ' (', makePct(stat$diploid   / stat$Total),    ')')
stat$Total     = paste0(stat$Total    , ' (', makePct(stat$Total     / stat$Total[1]), ')')
tfile          = paste0('1.infercnv.', name, '.predict.stat.xls')
write.table(cbind(Subcluster = rownames(stat), stat), tfile, sep = '\t', row.names = F, quote = F)
message(Sa('-->', timer(), 'saved', Pa(tfile), 'at', Pa(wdir), '<--'))
rm(stat)

# 2. Prepare statCNV #
limit  = quantile(ifcObj@expr.data, c(.01, .99))
pdata  = round(MinMax(ifcObj@expr.data, limit[1], limit[2]), 6)
Ref    = pdata[, unlist(ifcObj@reference_grouped_cell_indices)  , drop = F]
Obs    = pdata[, unlist(ifcObj@observation_grouped_cell_indices), drop = F]
rm(pdata)
## top annot:
tAnnot = ifcObj@gene_order[1]
#~~~~~~~~~~~~~
# ifcObj@gene_order:
#                 chr   start    stop
# ENSG00000228794   1  825138  859446
# ENSG00000188976   1  944203  959309
#~~~~~~~~~~~~~
## left annot:
lAnnot = obj@meta.data[plot_cols]
lAnnot = lAnnot[do.call(order, lAnnot[rev(names(lAnnot))]), , drop = F]
rm(ifcObj, obj)
# predict base-line #
num    = 1e+4
r_i    = sample(nrow(Ref), num, replace = T)
c_i    = sample(ncol(Ref), num, replace = T)
base   = as.numeric(names(sort(table( unlist(lapply(seq(num), function(n) 
  Ref[ r_i[n], c_i[n] ] ))), T)[1]))
rm(num, r_i, c_i)
# save
rfile  = paste0('0.infercnv.', name, '.statCNV.rds')
saveRDS(list(Obs = Obs - base, Ref = Ref - base, TopAnnot = tAnnot, LeftAnnot = lAnnot, CNV = CNV), rfile)
message(Sa('-->', timer(), 'saved', Pa(rfile), 'at', Pa(wdir), '<--'))

# done #
system('rm -f GD_inferCNV.undone ; touch GD_inferCNV.done')
message(Wa('==>', timer(), 'Done:', me, '<=='))

