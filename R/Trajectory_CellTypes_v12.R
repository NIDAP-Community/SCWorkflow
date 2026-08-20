#' Trajectory analysis with TSCAN
#'
#' Perform the CCBR single-cell RNA-seq trajectory workflow using TSCAN.
#'
#' @param Seurat_Object A Seurat object or path to an RDS file containing one.
#' @param MetaData A data.frame with per-cell metadata; must contain the clustering column.
#' @param Clusters_to_Use Character scalar naming the metadata column to use for clusters.
#' @param Custom_Gene_List Optional character vector or delimited string of genes to highlight.
#'
#' @return Invisibly returns `NULL`. Generates trajectory and pseudotime plots as side effects.
#' @examples
#' \donttest{
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   seurat_example <- Seurat::CreateSeuratObject(
#'     counts = matrix(rpois(200, lambda = 5), nrow = 20,
#'                     dimnames = list(paste0("g", 1:20), paste0("c", 1:10)))
#'   )
#'   seurat_example <- Seurat::AddMetaData(
#'     object = seurat_example,
#'     metadata = data.frame(cluster = sample(letters[1:2], 10, replace = TRUE),
#'                           row.names = colnames(seurat_example))
#'   )
#'   Trajectory_CellTypes(
#'     Seurat_Object = seurat_example,
#'     MetaData = Seurat::FetchData(seurat_example, vars = c("cluster")),
#'     Clusters_to_Use = "cluster"
#'   )
#' }
#' }
#' @export
Trajectory_CellTypes <- function(Seurat_Object,
                                 MetaData,
                                 Clusters_to_Use,
                                 Custom_Gene_List = NULL) {
  required_packages <- c(
    "Seurat", "scater", "scran", "scuttle", "TrajectoryUtils", "TSCAN",
    "ggplot2", "gridExtra", "mclust",
    "SingleCellExperiment", "SummarizedExperiment"
  )
  ensure_namespace(required_packages)

  seurat_object <- load_seurat_object(Seurat_Object)

  if (!is.data.frame(MetaData)) {
    stop("MetaData must be a data.frame.", call. = FALSE)
  }
  if (!is.character(Clusters_to_Use) || length(Clusters_to_Use) != 1L ||
      is.na(Clusters_to_Use)) {
    stop("Clusters_to_Use must be a non-missing character scalar.", call. = FALSE)
  }
  if (!Clusters_to_Use %in% colnames(MetaData)) {
    stop(sprintf("Column '%s' not found in MetaData.", Clusters_to_Use), call. = FALSE)
  }

  cache_original <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
  on.exit({
    if (is.na(cache_original)) {
      Sys.unsetenv("R_USER_CACHE_DIR")
    } else {
      Sys.setenv(R_USER_CACHE_DIR = cache_original)
    }
  }, add = TRUE)
  cache_dir <- tools::R_user_dir("Trajectory_CellTypes", which = "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(R_USER_CACHE_DIR = cache_dir)

  if (!is.null(rownames(MetaData)) && all(colnames(seurat_object) %in% rownames(MetaData))) {
    MetaData <- MetaData[colnames(seurat_object), , drop = FALSE]
  }

  user_clusters <- MetaData[[Clusters_to_Use]]
  if (length(user_clusters) != ncol(seurat_object)) {
    warning("MetaData rows do not match Seurat object cells; attempting to align by cell names.")
    shared_cells <- intersect(colnames(seurat_object), rownames(MetaData))
    if (length(shared_cells) == 0) {
      stop("Unable to align MetaData with Seurat object cells.", call. = FALSE)
    }
    MetaData <- MetaData[shared_cells, , drop = FALSE]
    MetaData <- MetaData[colnames(seurat_object), , drop = FALSE]
    user_clusters <- MetaData[[Clusters_to_Use]]
  }
  if (length(user_clusters) != ncol(seurat_object)) {
    stop("MetaData does not align with the Seurat object cells after attempted matching.", call. = FALSE)
  }

  label_column <- resolve_cluster_column(seurat_object, MetaData, Clusters_to_Use)
  if (!label_column %in% colnames(seurat_object@meta.data)) {
    seurat_object[[label_column]] <- user_clusters
  }
  cluster_labels <- seurat_object[[label_column, drop = TRUE]]
  Seurat::Idents(seurat_object) <- cluster_labels

  custom_gene_list <- normalise_gene_list(Custom_Gene_List)

  message("The following genes were used:")
  gene_universe <- rownames(seurat_object)
  genes_present <- intersect(custom_gene_list, gene_universe)
  if (length(genes_present) == 0) {
    message("None")
  } else {
    message(paste(genes_present, collapse = ", "))
  }

  message("The following genes were not found:")
  missing_genes <- setdiff(custom_gene_list, gene_universe)
  if (length(missing_genes) == 0) {
    message("None")
  } else {
    message(paste(missing_genes, collapse = ", "))
  }

  custom_gene_list <- genes_present

  sce <- Seurat::as.SingleCellExperiment(seurat_object)
  SummarizedExperiment::colData(sce)$label <- cluster_labels

  aggregated <- scuttle::aggregateAcrossCells(sce, ids = cluster_labels)
  centroids <- SingleCellExperiment::reducedDim(aggregated, "PCA")
  mst <- TSCAN::createClusterMST(centroids, clusters = NULL)

  line_data_tsne <- TSCAN::reportEdges(aggregated, mst = mst, clusters = NULL, use.dimred = "TSNE")
  line_data_tsne <- harmonise_edge_axes(line_data_tsne, "tSNE")
  plot_tsne1 <- scater::plotTSNE(sce, colour_by = "label") +
    ggplot2::geom_line(data = line_data_tsne, ggplot2::aes(x = tSNE_1, y = tSNE_2, group = edge))

  line_data_umap <- TSCAN::reportEdges(aggregated, mst = mst, clusters = NULL, use.dimred = "UMAP")
  line_data_umap <- harmonise_edge_axes(line_data_umap, "UMAP")
  plot_umap1 <- scater::plotUMAP(sce, colour_by = "label") +
    ggplot2::geom_line(data = line_data_umap, ggplot2::aes(x = UMAP_1, y = UMAP_2, group = edge))

  mapping <- TSCAN::mapCellsToEdges(sce, mst = mst, use.dimred = "PCA")
  ordered <- TSCAN::orderCells(mapping, mst)
  common_pseudo <- TrajectoryUtils::averagePseudotime(ordered)

  plot_tsne2 <- scater::plotTSNE(sce, colour_by = I(common_pseudo), text_by = "label", text_colour = "red") +
    ggplot2::geom_line(data = line_data_tsne, ggplot2::aes(x = tSNE_1, y = tSNE_2, group = edge)) +
    ggplot2::ggtitle("TSCAN-based MST")

  plot_umap2 <- scater::plotUMAP(sce, colour_by = I(common_pseudo), text_by = "label", text_colour = "red") +
    ggplot2::geom_line(data = line_data_umap, ggplot2::aes(x = UMAP_1, y = UMAP_2, group = edge)) +
    ggplot2::ggtitle("TSCAN-based MST")

  pseudo_all <- TSCAN::quickPseudotime(sce, use.dimred = "PCA")
  pseudo_mnn <- TSCAN::quickPseudotime(sce, use.dimred = "PCA", dist.method = "mnn")
  if (!is.null(pseudo_mnn$connected$TSNE)) {
    pseudo_mnn$connected$TSNE <- harmonise_edge_axes(pseudo_mnn$connected$TSNE, "tSNE")
  }
  if (!is.null(pseudo_mnn$connected$UMAP)) {
    pseudo_mnn$connected$UMAP <- harmonise_edge_axes(pseudo_mnn$connected$UMAP, "UMAP")
  }
  mnn_pseudo <- TrajectoryUtils::averagePseudotime(pseudo_mnn$ordering)

  plot_tsne3 <- scater::plotTSNE(sce, colour_by = I(mnn_pseudo), text_by = "label", text_colour = "red") +
    ggplot2::geom_line(data = pseudo_mnn$connected$TSNE, ggplot2::aes(x = tSNE_1, y = tSNE_2, group = edge)) +
    ggplot2::ggtitle("TSCAN with MNN distances-based MST")

  plot_umap3 <- scater::plotUMAP(sce, colour_by = I(mnn_pseudo), text_by = "label", text_colour = "red") +
    ggplot2::geom_line(data = pseudo_mnn$connected$UMAP, ggplot2::aes(x = UMAP_1, y = UMAP_2, group = edge)) +
    ggplot2::ggtitle("TSCAN with MNN distances-based MST")

  sce$Pseudotime <- TrajectoryUtils::pathStat(ordered)[, 1]
  pseudo_tests <- TSCAN::testPseudotime(sce, pseudotime = ordered[, 1])[[1]]
  sorted_tests <- pseudo_tests[order(pseudo_tests$p.value), ]

  downregulated <- sorted_tests[sorted_tests$logFC < 0, , drop = FALSE]
  downregulated$SYMBOL <- rownames(downregulated)
  top_down <- head(downregulated$SYMBOL, 10)
  plot_decrease <- scater::plotExpression(
    sce,
    features = top_down,
    x = "Pseudotime",
    colour_by = "label",
    ncol = 5
  ) + ggplot2::ggtitle("Expression of the top 10 genes that decrease with pseudotime")

  upregulated <- sorted_tests[sorted_tests$logFC > 0, , drop = FALSE]
  upregulated$SYMBOL <- rownames(upregulated)
  top_up <- head(upregulated$SYMBOL, 10)
  plot_increase <- scater::plotExpression(
    sce,
    features = top_up,
    x = "Pseudotime",
    colour_by = "label",
    ncol = 5
  ) + ggplot2::ggtitle("Expression of the top 10 genes that increase with pseudotime")

  on_first_path <- !is.na(sce$Pseudotime)

  print(Seurat::DimPlot(seurat_object))
  set.seed(1)
  plot(pseudo_all$mst)
  print(plot_tsne1)
  gridExtra::grid.arrange(plot_tsne2, plot_tsne3, ncol = 2)
  print(plot_umap1)
  gridExtra::grid.arrange(plot_umap2, plot_umap3, ncol = 2)
  plot(plot_decrease)
  plot(plot_increase)
  scater::plotHeatmap(
    sce[, on_first_path],
    order_columns_by = "Pseudotime",
    colour_columns_by = "label",
    features = head(top_up, 50),
    center = TRUE,
    main = "Expression of the top 50 genes that increase with pseudotime"
  )
  scater::plotHeatmap(
    sce[, on_first_path],
    order_columns_by = "Pseudotime",
    colour_columns_by = "label",
    features = head(top_down, 50),
    center = TRUE,
    main = "Expression of the top 50 genes that decrease with pseudotime"
  )

  if (length(custom_gene_list) > 0) {
    custom_scatter <- scater::plotExpression(
      sce,
      features = custom_gene_list,
      x = "Pseudotime",
      colour_by = "label",
      ncol = 5
    ) + ggplot2::ggtitle("Expression of custom genes across pseudotime")
    plot(custom_scatter)

    scater::plotHeatmap(
      sce[, on_first_path],
      order_columns_by = "Pseudotime",
      colour_columns_by = "label",
      features = custom_gene_list,
      center = TRUE,
      main = "Expression of custom genes across pseudotime"
    )
  }

  invisible(NULL)
}

utils::globalVariables(c(
  "label", "edge", "Pseudotime", "SYMBOL", "UMAP_1", "UMAP_2",
  "tSNE_1", "tSNE_2", "mst"
))

ensure_namespace <- function(pkgs) {
  missing <- vapply(pkgs, function(pkg) !requireNamespace(pkg, quietly = TRUE), logical(1))
  if (any(missing)) {
    stop(
      sprintf(
        "Missing required packages: %s",
        paste(pkgs[missing], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

load_seurat_object <- function(input) {
  if (inherits(input, "Seurat")) {
    return(input)
  }
  if (is.character(input) && length(input) == 1L && file.exists(input)) {
    return(readRDS(input))
  }
  if (requireNamespace("nidapFunctions", quietly = TRUE)) {
    path <- tryCatch(
      nidapFunctions::nidapGetPath(input, "seurat_object.rds"),
      error = function(e) NULL
    )
    if (!is.null(path) && file.exists(path)) {
      return(readRDS(path))
    }
  }
  stop(
    "Provide a Seurat object or a path to an RDS file containing one.",
    call. = FALSE
  )
}

resolve_cluster_column <- function(seurat_object, metadata, selected) {
  seurat_columns <- colnames(seurat_object@meta.data)
  if (selected %in% seurat_columns) {
    return(selected)
  }
  idx <- match(selected, colnames(metadata))
  if (!is.na(idx) && idx <= length(seurat_columns)) {
    candidate <- seurat_columns[idx]
    if (!is.na(candidate) && nzchar(candidate)) {
      return(candidate)
    }
  }
  normalised_selected <- make.names(selected)
  matches <- seurat_columns[make.names(seurat_columns) == normalised_selected]
  if (length(matches) == 1) {
    return(matches)
  }
  selected
}

normalise_gene_list <- function(genes) {
  if (is.null(genes) || length(genes) == 0) {
    return(character())
  }
  if (length(genes) == 1L) {
    genes <- unlist(strsplit(gsub(",", " ", genes), " "))
  }
  genes <- genes[genes != ""]
  unique(trimws(genes))
}

harmonise_edge_axes <- function(edge_df, prefix) {
  if (is.null(edge_df) || nrow(edge_df) == 0) {
    return(edge_df)
  }
  targets <- c(sprintf("%s_1", prefix), sprintf("%s_2", prefix))
  lower_targets <- tolower(targets)
  current_names <- names(edge_df)
  for (i in seq_along(targets)) {
    required <- targets[i]
    if (!required %in% current_names) {
      candidates <- current_names[tolower(current_names) == lower_targets[i]]
      if (length(candidates) == 1L) {
        names(edge_df)[names(edge_df) == candidates] <- required
      }
    }
  }
  edge_df
}
