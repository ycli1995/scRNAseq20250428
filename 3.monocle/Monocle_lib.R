############################################
############### subroutine #################

importCDS <- function (otherCDS, import_all = FALSE, assay.use = "RNA") {
		if (class(otherCDS)[1] == "seurat" || class(otherCDS)[1] == "Seurat") {
				requireNamespace("Seurat")
				if (grepl('^2', otherCDS@version)){
						data <- otherCDS@raw.data
						lowerDetectionLimit <- otherCDS@is.expr
				} else {
#				} else if (grepl('^3', otherCDS@version)) {
						if ( is.null(assay.use) ) assay.use <- DefaultAssay(otherCDS)
						data <- otherCDS@assays[[assay.use]]@counts
						if ( length(data) == 0 ) data <- otherCDS@assays[[assay.use]]@data
						lowerDetectionLimit <- 0
#				} else {
#						stop()
				}
				if (class(data) == "data.frame") {
						data <- as(as.matrix(data), "sparseMatrix")
				}
				pd <- tryCatch({
						pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
						pd
				}, error = function(e) {
						pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
						pd <- new("AnnotatedDataFrame", data = pData)
						message("This Seurat object doesn't provide any meta data")
						pd
				})
				if (length(setdiff(colnames(data), rownames(pd))) > 0) {
						data <- data[, rownames(pd)]
				}
				gene_short_name <- if ( exists("fdata", otherCDS@misc) ) FindFeaturesName(otherCDS, row.names(data)) else row.names(data)
				fData <- data.frame(gene_short_name = gene_short_name, 
						row.names = row.names(data))
				fd <- new("AnnotatedDataFrame", data = fData)
#				if (all(data == floor(data))) {  ## Error in asMethod(object) : Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105
				is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
				if (all(is.wholenumber(data[,head(1:ncol(data), 10000)]))) {
						expressionFamily <- negbinomial.size()
				}
				else if (any(data < 0)) {
#						expressionFamily <- gaussianff()
						expressionFamily <- uninormal()
				}
				else {
						expressionFamily <- tobit()
				}
				valid_data <- data[, row.names(pd)]
				monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
						lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
				if (import_all) {
						if ("Monocle" %in% names(otherCDS@misc)) {
								otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
								otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
								monocle_cds <- otherCDS@misc$Monocle
								mist_list <- otherCDS
						}
						else {
								mist_list <- otherCDS
						}
				}
				else {
						mist_list <- list()
				}
				if ("var.genes" %in% slotNames(otherCDS)) {
						var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
				}
				monocle_cds@auxClusteringData$seurat <- mist_list
		}
		else if (class(otherCDS)[1] == "SCESet") {
				requireNamespace("scater")
				message("Converting the exprs data in log scale back to original scale ...")
				data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
				fd <- otherCDS@featureData
				pd <- otherCDS@phenoData
				experimentData = otherCDS@experimentData
				if ("is.expr" %in% slotNames(otherCDS)) 
						lowerDetectionLimit <- otherCDS@is.expr
				else lowerDetectionLimit <- 1
						if (all(data == floor(data))) {
								expressionFamily <- negbinomial.size()
						}
						else if (any(data < 0)) {
								expressionFamily <- gaussianff()
						}
						else {
								expressionFamily <- tobit()
						}
				if (import_all) {
						mist_list <- otherCDS
				}
				else {
						mist_list <- list()
				}
				monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
						lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
				monocle_cds@auxOrderingData$scran <- mist_list
		}
		else {
				stop("the object type you want to export to is not supported yet")
		}
		return(monocle_cds)
}


Initial <- function( object, cells.use = NULL, assay.use = "RNA", feature_use = VariableFeatures(object) ) {
		if ( is.null(object) ){
				stop( "object must be defined." )
		}
		if ( is.null(cells.use) ) {
				select <- object
		}else {
				if (class(object) == "seurat") {
						select <- SubsetData( object, cells.use = cells.use, subset.raw = T )
				}else{
						select <- object[,cells.use]
				}
		}
		if ( is.null(assay.use) ) assay.use <- DefaultAssay(select)
		mnc_obj <- importCDS(select, assay.use = assay.use)
		fdata <- AddUnderscore(object@misc$fdata)
		new.name <- fdata$underscore
		names(new.name) <- fdata$dash
		index <- rownames(mnc_obj) %in% names(new.name)
		rownames(mnc_obj)[index] <- new.name[rownames(mnc_obj)[index]]
		mnc_obj <- estimateSizeFactors(mnc_obj)
		mnc_obj <- estimateDispersions(mnc_obj)
		if ( length(feature_use) == 0 ) {
			if ( mnc_obj@expressionFamily@vfamily %in% c('negbinomial', 'negbinomial.size') ){
				disp_table <- dispersionTable(mnc_obj)
				genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
			} else {
				stop(mnc_obj@expressionFamily@vfamily)
			}
		}else{
				feature_use <- new.name[feature_use]
				genes <- feature_use
		}
		mnc_obj <- setOrderingFilter(mnc_obj, ordering_genes = genes)
		
		return(mnc_obj)
}

select_ncenter <- function(mnc_obj, norm_method = c("log", "vstExprs", "none")) {
		if (ncol(mnc_obj) >= 100) {
				ncenter <- monocle:::cal_ncenter(ncol(mnc_obj))
		} else {
				ncenter <- ncol(mnc_obj) - 1
		}
		FM <- monocle:::normalize_expr_data(mnc_obj, norm_method, 1)
		xm <- Matrix::rowMeans(FM)
		xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
		FM <- FM[xsd > 0, ]
		FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
        FM <- FM[!is.na(row.names(FM)), ]
		FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ]
		X <- FM
		W <- DDRTree:::pca_projection_R(X %*% t(X), 2)
		Z <- t(W) %*% X
		ncenter <- seq(ncenter, 0)[!sapply(seq(ncenter, 0), function(i) any(duplicated(t(Z)[seq(1, ncol(Z), length.out = i), ])))][1]
		return(ncenter)
}


DoMonocle <- function( mnc_obj, ordering_genes = NULL, ... ) {
#				max_components = 2, maxIter = 20, ncenter = NULL, tol = 0.001, sigma = 0.001, lambda = NULL, param.gamma = 10 ) {
		if ( !is.null(ordering_genes) ) mnc_obj <- setOrderingFilter(mnc_obj, ordering_genes = ordering_genes)
		if ( !("negbinomial" == mnc_obj@expressionFamily@vfamily) || ("negbinomial.size" == mnc_obj@expressionFamily@vfamily) ) {
				norm_method <- "none"
		}else{
				norm_method <- c("log", "vstExprs", "none")
		}
		ncenter <- select_ncenter(mnc_obj, norm_method = norm_method)
		if ( length(ncenter) == 0 ){
				stop('unable to select proper ncenter, change ordering_genes may help.')
		}
		mnc_obj <- reduceDimension(mnc_obj, method = 'DDRTree', verbose = F, scaling = T, norm_method = norm_method, ncenter = ncenter, ... )
		mnc_obj <- orderCells(mnc_obj)
		return( mnc_obj )
}

findRC <- function(n) {
		row <- ceiling(sqrt(n))
		col <- ceiling(n/row)
		return(c(row, col))
}

Monocle_Plot <- function( mnc_obj, name, show_branch_points = FALSE, scale.size = 6 ) {
		facet.scale <- 0.7
		p1 <- plot_cell_trajectory(mnc_obj, color_by = "Samples",  show_branch_points = show_branch_points)
		rc <- findRC( length(levels(pData(mnc_obj)$Samples)) )
		p2 <- p1 + facet_wrap(~Samples,  nrow = rc[1])
		ggsave( paste0( name, "Sample.pdf" ),        p1, width = scale.size, height = scale.size, limitsize = FALSE )
		ggsave( paste0( name, "Sample.facet.pdf" ),  p2, width = scale.size * rc[2] * facet.scale, height = scale.size * rc[1] * facet.scale, limitsize = FALSE )

		p3 <- plot_cell_trajectory(mnc_obj, color_by = "Clusters", show_branch_points = show_branch_points)
		rc <- findRC( length(levels(pData(mnc_obj)$Clusters)) )
		p4 <- p3 + facet_wrap(~Clusters, nrow = rc[1])
		ggsave( paste0( name, "Cluster.pdf" ),       p3, width = scale.size, height = scale.size, limitsize = FALSE )
		ggsave( paste0( name, "Cluster.facet.pdf" ), p4, width = scale.size * rc[2] * facet.scale, height = scale.size * rc[1] * facet.scale, limitsize = FALSE )

		p5 <- plot_cell_trajectory(mnc_obj, color_by = "State",    show_branch_points = show_branch_points)
		rc <- findRC( length(levels(pData(mnc_obj)$State)) )
		p6 <- p5 + facet_wrap(~State,    nrow = rc[1])
		ggsave( paste0( name, "State.pdf" ),         p5, width = scale.size, height = scale.size, limitsize = FALSE )
		ggsave( paste0( name, "State.facet.pdf" ),   p6, width = scale.size * rc[2] * facet.scale, height = scale.size * rc[1] * facet.scale, limitsize = FALSE )

		p7 <- plot_cell_trajectory(mnc_obj, color_by = "Pseudotime", show_branch_points = show_branch_points )
		ggsave( paste0( name, "Pseudotime.pdf" ),    p7, width = scale.size, height = scale.size, limitsize = FALSE )

#		p6 <- plot_complex_cell_trajectory( mnc_obj, color_by = 'Cluster', show_branch_points = show_branch_points, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Sample, nrow = 1) + scale_size(range = c(0.2, 0.2))
#		ggsave( paste0( name, "tree_structure.pdf" ), width = 18, height = 6 )

#		saveRDS( mnc_obj, file = paste0( name, "obj.rds" ) )
}



select_sig_gene <- function( data, threshold = 1e-5, n = NULL, use.q = TRUE ) {
		if ( use.q ) {
				tmp <- subset( data[order(data$qval),], qval < threshold )
		}else{
				tmp <- subset( data[order(data$pval),], pval < threshold )
		}
		if ( is.null(n) ) {
				return( rownames(tmp) )
		}else {
				return( head(rownames(tmp), n) )
		}
}

FindBranchIndex <- function( cds ) {
		if (cds@dim_reduce_type == "DDRTree" ){
				reduced_dim_coords <- reducedDimK(cds)
		} else {
				stop("Error: unrecognized dimensionality reduction method.")
		}
		ica_space_df <- as.data.frame(Matrix::t(reduced_dim_coords))
		ica_space_df$pseudo_point <- rownames(ica_space_df)
		mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
		branch_point_df <- ica_space_df %>%
				slice(match(mst_branch_nodes, pseudo_point)) %>%
				mutate(branch_point_idx = seq_len(n())) 
		BranchIndex <- branch_point_df$branch_point_idx
		names(BranchIndex) <- branch_point_df$pseudo_point
		return(BranchIndex)
}

GetCellsInPath <- function(cds, branch_point = NULL, root_cell = NULL) {
		if ( is.null(root_cell) ) {
				root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
		}

		if(cds@dim_reduce_type == "DDRTree") {
				pr_graph_cell_proj_mst <- minSpanningTree(cds)
		}else {
				pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
		}

		closest_vertex <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex
		vertex2cell <- by(rownames(closest_vertex), paste0("Y_", closest_vertex), function(x) as.character(x))
		root_cell_Y <- paste0( "Y_", closest_vertex[root_cell,1] )
	
		path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_point, root_cell_Y)
		path_to_ancestor <- names(unlist(path_to_ancestor$vpath))

		mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_point]

		CellsInPath <- list()
		for (backbone_nei in neighbors( pr_graph_cell_proj_mst, branch_point )$name ) {
				descendents <- bfs(mst_no_branch_point, V(mst_no_branch_point)[backbone_nei], unreachable=FALSE)
				descendents <- descendents$order[!is.na(descendents$order)]
				descendents <- V(mst_no_branch_point)[descendents]$name
				if (root_cell_Y %in% descendents == FALSE){
						path_to_point <- unique( c(rev(path_to_ancestor), branch_point, descendents ) )
						if (cds@dim_reduce_type == "DDRTree"){
								cells_in_path <- unlist(vertex2cell[path_to_point], use.names = F)
						}else{
								cells_in_path <- path_to_point
						}
						cells_in_path <- intersect(cells_in_path, colnames(cds))
						CellsInPath[[backbone_nei]] <- cells_in_path
				}
		}
		return(CellsInPath)
}

.SummeriseCellsProperty <- function(cds, cells = NULL, col_name = 'State', collapse = TRUE) {
		if ( is.null(cells) ) {
				cells <- colnames(cds)
		}
		sub_p <- pData(cds)[cells, col_name, drop = F]
		p_name <- as.character(unique(sub_p[[col_name]]))
		if ( ! collapse ) {
				return(p_name)
		} else {
				name <- paste0("State ", paste(p_name, collapse = ","))
				return(name)
		}
}

FindBranchName <- function(cds, branch_point = NULL, root_cell = NULL) {
		cells_in_path <- GetCellsInPath(cds, branch_point, root_cell)
		BranchName <- lapply(cells_in_path, .SummeriseCellsProperty, cds = cds, col_name = "State")
		return(BranchName)
}


ChangeBranchPointState <- function(cds) {
		if(cds@dim_reduce_type == "DDRTree") {
				pr_graph_cell_proj_mst <- minSpanningTree(cds)
		}else {
				pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
		}

		closest_vertex      <- as.data.frame(cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex)
		closest_vertex$name <- paste0( "Y_", closest_vertex$V1 )
		root_cell      <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
		root_vertex    <- closest_vertex[root_cell, 2]

		for ( branch_point in cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points ) {
				path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_point, root_vertex)
				path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
				for ( i in 2:length(path_to_ancestor) ) {
						next_cells_to_root <- rownames(subset(closest_vertex, name == path_to_ancestor[i]))
						if ( length(next_cells_to_root) != 0 ) break
				}
				next_state_to_root <- unique(pData(cds)[next_cells_to_root, "State"])
				if ( length(next_state_to_root) == 1 ) {
						branch_cell <- rownames(subset(closest_vertex, name == branch_point))
						pData(cds)[branch_cell, "State"] <- next_state_to_root
				}else if ( length(next_state_to_root) == 0 ){
				}else{
						stop(next_state_to_root)
				}
		}

		return(cds)
}


############### __subroutine__ #################
################################################





MonocleObject <- function(parameter = list(),
				cls_col_name = parameter$cluster_col, sample_col_name = parameter$sample_col,
				mnc_obj = parameter$mnc_obj, obj_file = parameter$obj_file,
				cell_use = parameter$cell_use, #feature_use = parameter$feature_use,
				sample_use = parameter$sample_use, samele_exclude = parameter$samele_exclude,
				cluster_use = parameter$cluster_use, cluster_exclude = parameter$cluster_exclude,
				rename_samples = parameter$samples, rename_clusters = parameter$rename_clusters,
				cellSample = parameter$cellSample, assay.use = parameter$assay.use
				){
		if ( is.null(cls_col_name) ) cls_col_name = "seurat_clusters"
		if ( is.null(sample_col_name) ) sample_col_name = "orig.ident"
		if ( ! is.null(mnc_obj) && file.exists(mnc_obj) ) {
				mnc_obj <- Load(mnc_obj)
		} else {
				library(Seurat)
				### read obj 
				obj <- Load(obj_file)

				if ( ! is.null(cell_use) ) {
						if ( ! is.list(cell_use) ) {
								cells <- readLines(cell_use)
								cells <- intersect(cells, Cells(obj))
								if ( length(cells) == 0 ) {
										stop(cell_use, " doesnot intersect with Cells(obj)")
								}
								obj <- obj[, cells]
						} else {
						cellset <- ReadCellSet(cell_use)

						cellinfo <- read.table(cellSample, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
#						orig_name <- cellinfo$orig_name
#						names(orig_name) <- rownames(cellinfo)

						if ( all(cellset$cells %in% Cells(obj)) ) {
								message("All select cells are in the object. No need to RemakeObject.")
								obj <- CellSetFilter(obj, cell_use, col.name = cls_col_name)
#								if ( ! identical(as.vector(orig_name[Cells(obj)]), Cells(obj)) ) {
#										obj <- RenameCells(obj, new.names = as.vector(orig_name[Cells(obj)]))
#								}       
								obj@meta.data[[sample_col_name]] <- cellinfo[Cells(obj), "sample"]
						} else {
								obj <- RemakeObject(obj)
								now_name <- rownames(cellinfo)
								names(now_name) <- cellinfo$orig_name
								tmp <- Cells(obj)[!Cells(obj) %in% names(now_name)]
								names(tmp) <- tmp
								now_name <- c(now_name, tmp)
								obj <- RenameCells(obj, new.names = as.vector(now_name[Cells(obj)]))
								obj <- CellSetFilter(obj, cell_use, col.name = cls_col_name)
						}
						}
				}
				if ( ! is.null(rename_samples) ) {
						new_sample <- NULL
						new_group  <- NULL
						for ( i in names(rename_samples) ) {
								new_sample[[i]] <- rename_samples[[i]][1]
								new_group[[rename_samples[[i]][1]]] <- rename_samples[[i]][2]
						}
						obj <- RenameGrouping(obj, new_sample = new_sample, new_group = new_group)
				}
				if ( ! is.null(rename_clusters) ) {
						new_cluster <- unlist.rev(rename_clusters)
						index <- levels(obj@meta.data[[cls_col_name]]) %in% names(new_cluster)
						levels(obj@meta.data[[cls_col_name]])[index] <- new_cluster[levels(obj@meta.data[[cls_col_name]])[index]]
						Idents(obj) <- cls_col_name
						obj@misc$color.cluster <- SetColor(Idents(obj))
				}


				obj@meta.data$Samples  <- obj@meta.data[[sample_col_name]]
				if ( ! is.factor(obj@meta.data$Samples) ) obj@meta.data$Samples <- as.factor(obj@meta.data$Samples)
				obj@meta.data$Clusters <- obj@meta.data[[cls_col_name]]
				if ( ! is.factor(obj@meta.data$Clusters) ) obj@meta.data$Clusters <- as.factor(obj@meta.data$Clusters)


				### filtering
				sample_use  <- if ( ! is.null(sample_use)  ) sample_use else levels(obj@meta.data$Samples)
				sample_use  <- setdiff( sample_use, samele_exclude )
				cluster_use <- if ( ! is.null(cluster_use) ) cluster_use else levels(obj@meta.data$Clusters)
				cluster_use <- setdiff( cluster_use, cluster_exclude )
				cells.use   <- rownames(obj@meta.data)[obj@meta.data$Samples %in% sample_use & obj@meta.data$Clusters %in% cluster_use]

				### covert seurat to monocle
				mnc_obj <- Initial( obj, cells.use, assay.use)

		}
		return(mnc_obj)
}


RunMonocle <- function(mnc_obj, parameter = list(), feature_use = parameter$feature_use, ...){
				genes.use <- if ( ! is.null(feature_use) && file.exists(feature_use) ) readLines(feature_use) else NULL
				mnc_obj <- DoMonocle( mnc_obj, ordering_genes = genes.use, ... )
				pData(mnc_obj) <- droplevels(pData(mnc_obj))
				return(mnc_obj)
}

RenameState <- function(mnc_obj, parameter = list(), rename_state = parameter$rename_state, root_state = parameter$root_state, reverse = NULL, force = FALSE){
		if ( !is.null(root_state) ) {
				if ( force || pData(mnc_obj)$State[which.min(pData(mnc_obj)$Pseudotime)] != root_state ) {
						if ( ! is.null(reverse) && reverse ) root_state <- NULL
						mnc_obj <- orderCells( mnc_obj, root_state = root_state, reverse = reverse )
				}
		}
		if ( ! is.null(rename_state) ){
				if ( length(rename_state) < length(levels(pData(mnc_obj)$State)) ) stop("length of <rename_state> is less than State's levels.") 
				levels(pData(mnc_obj)$State) <- rename_state
		}
		return(mnc_obj)
}


StatTrajectory <- function(mnc_obj, name = "Trajectory.", scale.size = 6, show_branch_points = FALSE){
		Monocle_Plot( mnc_obj, name = name, scale.size = scale.size, show_branch_points = show_branch_points )

		Trajectory.data <- pData(mnc_obj) %>%
				tibble::rownames_to_column(var = "cells.name") %>%
				select(cells.name, Samples, Clusters, Pseudotime, State)
		write.table( Trajectory.data, file = paste0(name, "data.xls"), quote = F, sep = "\t", row.names = F )
}


Diff.State <- function(mnc_obj, thres = 1e-7, cores = 1, use.q = TRUE, diff_state_res = NULL){
		if ( is.null(diff_state_res) ) {
			diff_state_res <- differentialGeneTest(mnc_obj, fullModelFormulaStr = "~State", cores = cores)
		}
		saveRDS( diff_state_res, file = "diff_state.rds" )

		thres <- as.numeric(thres)
		sig_gene_names <- select_sig_gene(diff_state_res, thres, use.q = use.q)

		if ( length(sig_gene_names) > 0 ) {
		p1 <- plot_genes_jitter(mnc_obj[head(sig_gene_names, 10),], grouping = "State", color_by = "State", ncol = 5)
		ggsave( "Diff.state.pdf", p1, width = 12, height = 5, limitsize = FALSE )

		## for online report
		nstates <- length(levels(pData(mnc_obj)$State))
		nstates <- min(nstates, length(head(sig_gene_names,1000)))
		p3 <- plot_pseudotime_heatmap(mnc_obj[head(sig_gene_names,1000),], num_clusters = nstates, cores = cores, show_rownames = T, return_heatmap = T)
		ggsave( p3[["ph_res"]], file = "Diff.pseudotime_heatmap.pdf", width = 7, height = min(30, max(7, length(sig_gene_names) * 0.1)), limitsize = FALSE ) ## file name !
		write.table(sig_gene_names, file = "Diff.state_heatmap.list.xls", quote = F, sep = "\t", row.names = F, col.names = F)


		tmp <- by( rownames(pData(mnc_obj)), pData(mnc_obj)$State, function(x) Matrix::rowMeans(Biobase::exprs(mnc_obj)[sig_gene_names, x, drop = F]) )
		tmp <- do.call(cbind.data.frame, tmp)
		colnames(tmp) <- paste0( "State", colnames(tmp) )
		diff_state_sig <- diff_state_res[sig_gene_names,]
		diff_state_sig <- data.frame( GeneID = rownames(diff_state_sig), tmp, P_value = diff_state_sig$pval, FDR = diff_state_sig$qval, row.names = rownames(diff_state_sig) )
		write.table( diff_state_sig, file = "Diff.state.xls",  quote = F, sep = "\t", row.names = F )

		invisible(diff_state_res[sig_gene_names,])
		}
}

Diff.Pseudotime <- function(mnc_obj, thres = 1e-7, cores = 1, use.q = TRUE, diff_Pseudotime_res = NULL){
		if ( is.null(diff_Pseudotime_res) ) {
			diff_Pseudotime_res  <- differentialGeneTest(mnc_obj, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = cores)
			diff_Pseudotime_res  <- diff_Pseudotime_res[rownames(mnc_obj),]
			diff_Pseudotime_res  <- diff_Pseudotime_res[! is.na(diff_Pseudotime_res[[1]]),]
		}
		saveRDS( diff_Pseudotime_res, file = "diff_Pseudotime.rds" )

		thres <- as.numeric(thres)
		sig_gene_names <- select_sig_gene(diff_Pseudotime_res, thres = thres, use.q = use.q)

		if ( length(sig_gene_names) > 0 ) {
		p2 <- plot_genes_in_pseudotime(mnc_obj[head(sig_gene_names,10),], ncol = 5, color_by = "Samples")
		ggsave("Diff.genes_in_pseudotime.pdf", p2, width = 12, height = 5, limitsize = FALSE )

		nstates <- length(levels(pData(mnc_obj)$State))
		p3 <- plot_pseudotime_heatmap(mnc_obj[sig_gene_names,], num_clusters = nstates, cores = cores, show_rownames = T, return_heatmap = T)
		ggsave( p3[["ph_res"]], file = "Diff.pseudotime_heatmap.pdf", width = 7, height = min(30, max(7, length(sig_gene_names) * 0.1)), limitsize = FALSE )
		write.table(sig_gene_names, file = "Diff.pseudotime_heatmap.list.xls", quote = F, sep = "\t", row.names = F, col.names = F)

		cluster <- cutree(p3[["ph_res"]]$tree_row, nstates)[p3[["ph_res"]]$tree_row$order]
		diff_Pseudotime_sig <- diff_Pseudotime_res[intersect(sig_gene_names, names(cluster)),]
		diff_Pseudotime_sig <- data.frame( GeneID = rownames(diff_Pseudotime_sig), Gene_cluster = cluster[rownames(diff_Pseudotime_sig)], P_value = diff_Pseudotime_sig$pval, FDR = diff_Pseudotime_sig$qval, row.names = rownames(diff_Pseudotime_sig) )
		write.table( diff_Pseudotime_sig, file = "Diff.pseudotime.xls", quote = F, sep = "\t", row.names = F )

		invisible(diff_Pseudotime_res[sig_gene_names,])
		}
}


Diff.Branch <- function(mnc_obj, thres = 1e-7, cores = 1, use.q = TRUE, BEAM_res = NULL) {
		branch_point <- mnc_obj@auxOrderingData[[mnc_obj@dim_reduce_type]]$branch_points
		branch_point_index <- FindBranchIndex( mnc_obj )

		if ( is.null(BEAM_res) ) BEAM_res <- list()
		total_BEAM_sig <- list()
		for ( i in seq(branch_point) ) {
				cells_in_path <- GetCellsInPath(mnc_obj, branch_point[i])
				common_cell <- Reduce(intersect, cells_in_path)
				if ( length(cells_in_path) < 2 || any(sapply(cells_in_path, function(i) length(setdiff(i, common_cell))) == 0) ) {
						next
				}
				branch_name <- FindBranchName( ChangeBranchPointState(mnc_obj), branch_point[i] )
				if ( ! i %in% names(BEAM_res) ) {
						beam_res <- BEAM(mnc_obj, branch_point = i, cores = cores, branch_labels = branch_name )
						BEAM_res[[as.character(i)]] <- beam_res
				} else {
						beam_res <- BEAM_res[[as.character(i)]]
				}
				saveRDS( beam_res, file = paste0( "BEAM.", i, ".rds") )

				thres <- as.numeric(thres)
				sig_gene_names <- select_sig_gene(beam_res, thres = thres, use.q = use.q)

				if ( length(sig_gene_names) > 0 ) {
				p4 <- plot_genes_branched_pseudotime(mnc_obj[head(sig_gene_names,10),], branch_point = i, color_by = "Samples", ncol = 5, branch_labels = branch_name)
				ggsave( p4, file = paste0( "Branch.", branch_point_index[branch_point[i]], ".genes_pseudotime.pdf" ), width = 12, height = 5, limitsize = FALSE )

				num_clusters <- min(length(sig_gene_names), 5)
				p5 <- plot_genes_branched_heatmap(mnc_obj[sig_gene_names,], branch_point = i, num_clusters = num_clusters, cores = cores, show_rownames = T, return_heatmap=T, branch_labels = branch_name )
				ggsave( p5[["ph_res"]], file = paste0( "Branch.", branch_point_index[branch_point[i]], ".genes_heatmap.pdf" ), width = 7, height = min(30, max(7, length(sig_gene_names) * 0.1)), limitsize = FALSE )
				write.table(sig_gene_names, file = paste0( "Branch.", branch_point_index[branch_point[i]], ".genes_heatmap.list.xls"), quote = F, sep = "\t", row.names = F, col.names = F)

				if ( ! is.null(p5$annotation_row) ) {
						BEAM_sig <- beam_res[beam_res$gene_short_name %in% rownames(p5$annotation_row),,drop = FALSE]
						Gene_cluster <- p5$annotation_row[as.character(BEAM_sig$gene_short_name), "Cluster"]
				} else {
						BEAM_sig <- beam_res[sig_gene_names,,drop = FALSE]
						Gene_cluster <- '-'
				}
				BEAM_sig <- data.frame(GeneID = rownames(BEAM_sig), Gene_cluster = Gene_cluster, P_value = BEAM_sig$pval, FDR = BEAM_sig$qval, row.names = rownames(BEAM_sig))
				write.table( BEAM_sig, file = paste0( "Branch.", branch_point_index[branch_point[i]], ".depended_gene.xls" ), quote = F, sep = "\t", row.names = F )

				total_BEAM_sig[[i]] <- data.frame(branch_node = branch_point_index[branch_point[i]],
				                                  branch      = paste(branch_name, collapse = " -vs- "),
												  BEAM_sig)
				}
		}
		saveRDS(BEAM_res, file = "BEAM.rds")

		total_BEAM_sig <- do.call(rbind, total_BEAM_sig)
		write.table( total_BEAM_sig, file = "Branch.depended_gene.xls", quote = F, sep = "\t", row.names = F )

		invisible(total_BEAM_sig)
}


GetTrajectoryData <- function (cds, x = 1, y = 2, theta = 0, is.return = FALSE) {
#    requireNamespace("igraph")
    sample_state <- pData(cds)$State
    lib_info_with_pseudo <- pData(cds)
    if (is.null(cds@dim_reduce_type)) {
        stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
    }
    if (cds@dim_reduce_type == "ICA") {
        reduced_dim_coords <- reducedDimS(cds)
    } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
        reduced_dim_coords <- reducedDimK(cds)
    } else {
        stop("Error: unrecognized dimensionality reduction method.")
    }
    ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
        select(prin_graph_dim_1 := !!x, prin_graph_dim_2 := !!y) %>% 
        mutate(sample_name = rownames(.), sample_state = rownames(.))
    dp_mst <- minSpanningTree(cds)
    if (is.null(dp_mst)) {
        stop("You must first call orderCells() before using this function")
    }
    edge_df <- dp_mst %>% igraph::as_data_frame() %>% select(source = from, 
        target = to) %>% left_join(ica_space_df %>% select(source = sample_name, 
        source_prin_graph_dim_1 = prin_graph_dim_1, source_prin_graph_dim_2 = prin_graph_dim_2), 
        by = "source") %>% left_join(ica_space_df %>% select(target = sample_name, 
        target_prin_graph_dim_1 = prin_graph_dim_1, target_prin_graph_dim_2 = prin_graph_dim_2), 
        by = "target")
    data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
        select(data_dim_1 := !!x, data_dim_2 := !!y) %>% tibble::rownames_to_column("sample_name") %>% 
        mutate(sample_state) %>% left_join(lib_info_with_pseudo %>% 
        tibble::rownames_to_column("sample_name"), by = "sample_name")

if ( theta!= 0 ){
    return_rotation_mat <- function(theta) {
        theta <- theta/180 * pi
        matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
            nrow = 2)
    }
    rot_mat <- return_rotation_mat(theta)
    cn1 <- c("data_dim_1", "data_dim_2")
    cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
    cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
    data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
    edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
    edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
}

	branch_point_index <- FindBranchIndex(cds)
	edge_df <- edge_df %>% mutate(
			source_branch_point_index = if_else(is.na(branch_point_index[source]), "NULL", as.character(branch_point_index[source])),
			target_branch_point_index = if_else(is.na(branch_point_index[target]), "NULL", as.character(branch_point_index[target])))
	data_df <- data_df %>% select(sample_name, data_dim_1, data_dim_2, Samples, Clusters, Pseudotime, State, Size_Factor)

	if ( is.return ) {
			return(list(cells = data_df, bone = edge_df))
	}else{
			write.table(data_df, file = "Trajectory.cell.data.xls", sep = "\t", quote = F, row.names = F)
			write.table(edge_df, file = "Trajectory.bone.data.xls", sep = "\t", quote = F, row.names = F)
	}
}


GetGenesInPseudotime <- function (cds_subset, min_expr = NULL, trend_formula = "~ sm.ns(Pseudotime, df=3)", cores = 1, is.return = FALSE, outfile = NULL ) {
	if ( is.null(outfile) ) outfile <- "Trajectory.genes.pseudotime.data.xls" 
	cds_exprs <- exprs(cds_subset)
	if ( cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
			cds_exprs_rel <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
	}else{
			cds_exprs_rel <- cds_exprs
	}
	cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
	cds_exprs_rel <- reshape2::melt(round(as.matrix(cds_exprs_rel)))
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
	colnames(cds_exprs_rel) <- c("f_id", "Cell", "expression_relative")
	cds_exprs <- merge(cds_exprs, cds_exprs_rel)
    cds_exprs <- merge(cds_exprs, fData(cds_subset), by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, pData(cds_subset), by.x = "Cell", by.y = "row.names")

    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
	expectation <- reshape2::melt(model_expectation, varnames = c("f_id", "Cell"), value.name = "expectation")
    cds_exprs <- merge(cds_exprs, expectation)

	if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr

	cds_exprs <- cds_exprs %>% select(Cell, GeneID = f_id, GeneName = gene_short_name, expression, expression_relative, expectation, Pseudotime, State, Samples, Clusters)

	if ( is.return ) {
			return(cds_exprs)
	}else{
			write.table(cds_exprs, file = outfile, sep = "\t", quote = F, row.names = F)
	}
}



GetGenesBranchedPseudotime <- function (cds, branch_point = 1, branch_labels = NULL, method = "fitting", min_expr = NULL, trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch", reducedModelFormulaStr = NULL, is.return = FALSE, cores = 1, outfile = NULL, ...) {
	if ( is.null(outfile) ) outfile <- "Trajectory.genes.branch.data.xls"
	branch_states <- NULL # don't know what it used for yet
#    if (is.null(reducedModelFormulaStr) == FALSE) {
#        pval_df <- branchTest(cds, branch_states = branch_states, branch_point = branch_point, fullModelFormulaStr = trend_formula, reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", ...)
#        fData(cds)[, "pval"] <- pval_df[row.names(cds), "pval"]
#    }
    if ("Branch" %in% all.vars(terms(as.formula(trend_formula)))) {
        cds_subset <- buildBranchCellDataSet(cds = cds, branch_states = branch_states, branch_point = branch_point, branch_labels = branch_labels, progenitor_method = "duplicate", ...)
    } else {
        cds_subset <- cds
        pData(cds_subset)$Branch <- pData(cds_subset)$State
    }

	cds_exprs <- exprs(cds_subset)
    if ( cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size") ) {
        cds_exprs_rel <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    } else {
		cds_exprs_rel <- cds_exprs
    }
	cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
	cds_exprs_rel <- reshape2::melt(round(as.matrix(cds_exprs_rel)))
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
	colnames(cds_exprs_rel) <- c("f_id", "Cell", "expression_relative")
	cds_exprs <- merge(cds_exprs, cds_exprs_rel)
    cds_exprs <- merge(cds_exprs, fData(cds_subset), by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, pData(cds_subset), by.x = "Cell", by.y = "row.names")
    cds_exprs$Branch <- as.factor(cds_exprs$Branch)

    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)
    full_model_expectation <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, relative_expr = T, new_data = new_data)
    colnames(full_model_expectation) <- colnames(cds_subset)
	expectation <- reshape2::melt(full_model_expectation, varnames = c("f_id", "Cell"), value.name = "full_model_expectation")
	cds_exprs <- merge(cds_exprs, expectation)
#    if (!is.null(reducedModelFormulaStr)) {
#        reduced_model_expectation <- genSmoothCurves(cds_subset, cores = cores, trend_formula = reducedModelFormulaStr, relative_expr = T, new_data = new_data)
#        colnames(reduced_model_expectation) <- colnames(cds_subset)
#        cds_exprs$reduced_model_expectation <- apply(cds_exprs, 1, function(x) reduced_model_expectation[x[2], x[1]])
#    }
    if (method == "loess") cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit

    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
    cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
#    if (!is.null(reducedModelFormulaStr)) {
#        cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
#        cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
#    }
    cds_exprs$State <- as.factor(cds_exprs$State)
    cds_exprs$Branch <- as.factor(cds_exprs$Branch)

	cds_exprs <- cds_exprs %>% select(Cell, GeneID = f_id, GeneName = gene_short_name, expression, expression_relative, expectation = full_model_expectation, Pseudotime, State, Samples, Clusters, original_cell_id, Branch)

	if ( is.return ){
			return(cds_exprs)
	}else{
			write.table(cds_exprs, file = outfile, sep = "\t", quote = F, row.names = F)
	}

}


PlotPseudotimeHeatmap <- function(mnc_obj, features, parameter, cores = 1, outfile = NULL,
				position = parameter$legend$position,
				title = parameter$labs$title,
				color = unlist(parameter$legend$color),
				cluster_rows = parameter$hclust$row,
				scale.range = parameter$data$scale.range,
				xlab = parameter$labs$xlab, ylab = parameter$labs$ylab,
				show_rownames = parameter$axis$y.text$show, show_colnames = parameter$axis$x.text$show,
				fontsize = parameter$font$size, rowname_type = parameter$data$rowname_type,
				...){
		is.legend <- if(position == "none") FALSE else TRUE
		title <- if (is.null(title)) NA else title
		nstates <- length(levels(pData(mnc_obj)$State))

		if ( length(unlist(color)) < 3 ) {
				hmcols <- monocle:::blue2green2red(100)
		}else{
				hmcols <- colorRampPalette(color)(100)
		}
		ph <- plot_pseudotime_heatmap(mnc_obj[features,], num_clusters = nstates, cores = cores, show_rownames = T, return_heatmap = T,
						use_gene_short_name = FALSE,
						cluster_rows = cluster_rows, hmcols = hmcols, scale_max = scale.range, scale_min = scale.range * -1)
		height <- grid::convertHeight(sum(ph$ph_res$gtable$heights), "inches", valueOnly = T)
		width  <- grid::convertWidth(sum(ph$ph_res$gtable$widths), "inches", valueOnly = T)
		while(!is.null(dev.list())) { 
				dev.off()
		}

		if ( !is.null(outfile) ) {
				pdf(file = outfile,
					height = height, width = width)
		}
		if ( ! is.null(xlab) || ! is.null(ylab) ) {
				require(grid)
				setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
		}
		if ( ! is.null(ph$annotation_row) ) colnames(ph$annotation_row) <- "Gene_cluster"
		mat <- ph$heatmap_matrix[,]
		if ( is.null(rowname_type) || rowname_type == "id" ) {
		} else if ( rowname_type == "name" ) {
				name <- fData(mnc_obj)[rownames(mat), "gene_short_name"]
				rownames(mat) <- name
				if ( ! is.null(ph$annotation_row) ) rownames(ph$annotation_row) <- name
		} else if ( rowname_type == "id_name" ) {
				name <- fData(mnc_obj)[rownames(mat), "gene_short_name"]
				name <- paste0(rownames(mat), "(", name, ")")
				rownames(mat) <- name
				if ( ! is.null(ph$annotation_row) ) rownames(ph$annotation_row) <- name
		}
		pheatmap::pheatmap(mat,
				useRaster = T, border_color = NA, 
				color = ph$hmcols, breaks = ph$bks,
				clustering_method = ph$clustering_method,
				cutree_rows = ph$num_clusters, 
				cluster_rows = cluster_rows, cluster_cols = FALSE,
				clustering_distance_rows = ph$row_dist,
				treeheight_row = 20, 
				legend = is.legend,
				annotation_row = ph$annotation_row, annotation_col = ph$annotation_col,
				annotation_names_row = FALSE, annotation_names_col = FALSE,
				show_rownames = show_rownames, show_colnames = show_colnames,
				main = title,
				fontsize = fontsize,
				labels_row = NULL,
				...)
		if ( ! is.null(xlab) || ! is.null(ylab) ) {
				setHook("grid.newpage", NULL, "replace")
				grid.text(xlab, y=-0.05, gp=gpar(fontsize=fontsize))
				grid.text(ylab, x=-0.05, rot=90, gp=gpar(fontsize=fontsize))
		}

		if ( !is.null(outfile) ) {
				dev.off()
		}
}

PlotBranchHeatmap <- function(mnc_obj, features, parameter, cores = 1, outfile = NULL, branch_point = NULL,
				position = parameter$legend$position,
				title = parameter$labs$title,
				color = unlist(parameter$legend$color),
				scale.range = parameter$data$scale.range,
				cluster_rows = parameter$hclust$row,
				xlab = parameter$labs$xlab, ylab = parameter$labs$ylab,
				show_rownames = parameter$axis$y.text$show, show_colnames = parameter$axis$x.text$show,
				fontsize = parameter$font$size, rowname_type = parameter$data$rowname_type,
				...){
		is.legend <- if(position == "none") FALSE else TRUE
		title <- if (is.null(title)) NA else title

		branch_point_index <- FindBranchIndex( mnc_obj )
		branch_point_name <- names(branch_point_index)[branch_point_index == branch_point]
		branch_point <- branch_point_index[branch_point_index == branch_point]
		branch_name <- FindBranchName(ChangeBranchPointState(mnc_obj), branch_point_name)

		num_clusters <- 5
		if ( length(color) < 3 ) {
				hmcols <- monocle:::blue2green2red(100)
		}else{
				hmcols <- colorRampPalette(color)(100)
		}
		ph <- plot_genes_branched_heatmap(mnc_obj[features,], branch_point = branch_point,
						num_clusters = num_clusters, cores = cores,
						show_rownames = T, return_heatmap=T,
						branch_labels = branch_name, cluster_rows = cluster_rows, hmcols = hmcols,
						use_gene_short_name = FALSE,
						scale_max = scale.range,
						scale_min = scale.range * -1)
		height <- grid::convertHeight(sum(ph$ph_res$gtable$heights), "inches", valueOnly = T)
		width  <- grid::convertWidth(sum(ph$ph_res$gtable$widths), "inches", valueOnly = T)
		while(!is.null(dev.list())) { 
				dev.off()
		}

		if ( ! is.null(outfile) ) {
				pdf(file = outfile,
					height = height, width = width)
		}
		if ( ! is.null(xlab) || ! is.null(ylab) ) {
				require(grid)
				setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
		}
		if ( ! is.null(ph$annotation_row) ) colnames(ph$annotation_row) <- "Gene_cluster"
		mat <- ph$heatmap_matrix[,]
		if ( is.null(rowname_type) || rowname_type == "id" ) {
		} else if ( rowname_type == "name" ) {
				name <- fData(mnc_obj)[rownames(mat), "gene_short_name"]
				rownames(mat) <- name
				if ( ! is.null(ph$annotation_row) ) rownames(ph$annotation_row) <- name
		} else if ( rowname_type == "id_name" ) {
				name <- fData(mnc_obj)[rownames(mat), "gene_short_name"]
				name <- paste0(rownames(mat), "(", name, ")")
				rownames(mat) <- name
				if ( ! is.null(ph$annotation_row) ) rownames(ph$annotation_row) <- name
		}
		pheatmap::pheatmap(mat,
				useRaster = T, border_color = NA, 
				color = ph$hmcols, breaks = ph$bks,
				clustering_method = "ward.D2",
				cutree_rows = num_clusters,
				cluster_rows = cluster_rows, cluster_cols = FALSE,
				clustering_distance_rows = ph$row_dist,
				treeheight_row = 20, 
				legend = is.legend,
				annotation_row = ph$annotation_row, annotation_col = ph$annotation_col,
				annotation_colors = ph$annotation_colors,
				annotation_names_row = FALSE, annotation_names_col = FALSE,
				show_rownames = show_rownames, show_colnames = show_colnames,
				main = title,
				fontsize = fontsize,
				labels_row = NULL,
				gaps_col = ph$col_gap_ind,
#				filename = outfile, width = 7, height = NA,
				...)
		if ( ! is.null(xlab) || ! is.null(ylab) ) {
				setHook("grid.newpage", NULL, "replace")
				grid.text(xlab, y=-0.05, gp=gpar(fontsize=fontsize))
				grid.text(ylab, x=-0.05, rot=90, gp=gpar(fontsize=fontsize))
		}

		if ( !is.null(outfile) ) {
				dev.off()
		}
}

CDS_avg <- function(cds, group.by = "State", outfile = "AllGene.avg_exp.xls"){
		tmp <- by( rownames(pData(cds)), pData(cds)[[group.by]], function(x) Matrix::rowMeans(Biobase::exprs(cds)[, x, drop = F]) )
		tmp <- do.call(cbind.data.frame, tmp)
		colnames(tmp) <- paste0( group.by, colnames(tmp) )
		tmp <- data.frame( GeneID = rownames(tmp), GeneName = fData(cds)$gene_short_name, tmp)
		WriteTable(tmp, file = outfile )
}



