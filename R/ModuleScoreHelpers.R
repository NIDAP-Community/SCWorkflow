#' @title Helpers for ModuleScore Shiny app
#' @description Precompute module scores per celltype and build plots from cached data.
#' @name modscore-imports
#' @keywords internal
#' @importFrom dplyr mutate group_by summarise arrange select
#' @importFrom ggplot2 ggplot aes theme_bw theme element_blank 
#' @importFrom ggplot2 geom_point scale_color_gradientn guides guide_legend 
#' @importFrom ggplot2 xlab ylab element_text geom_violin theme_classic geom_hline 
#' @importFrom ggplot2 scale_y_continuous geom_line geom_segment scale_y_log10 scale_x_continuous
#' @importFrom gridExtra arrangeGrob 
#' @importFrom grid textGrob gpar
NULL

#' @export
compute_modscore_data <- function(object, marker.list, use.columns,
                                  reduction = c("tsne","umap","pca"),
                                  nbins = 10,
                                  group.var = "orig.ident") {
  reduction <- match.arg(reduction)

  if (!"Barcode" %in% colnames(object@meta.data)) {
    object@meta.data$Barcode <- rownames(object@meta.data)
  }
  colnames(object@meta.data) <- gsub("orig_ident", "orig.ident", colnames(object@meta.data))

  # Prepare identities and embeddings once
  idents <- Seurat::Idents(object)
  emb <- Seurat::Embeddings(object, reduction)
  if (ncol(emb) < 2) stop("Selected reduction has fewer than 2 dimensions.")

  # Build per-celltype data
  res <- list()
  for (celltype_name in use.columns) {
    genes <- marker.list[[celltype_name]]
    if (is.null(genes)) next

    object <- Seurat::AddModuleScore(object, list(genes), name = celltype_name, nbin = nbins, assay = "SCT")
    m <- paste0(celltype_name, "1")
    object@meta.data[[m]] <- scales::rescale(object@meta.data[[m]], to = c(0, 1))

    clusid <- as.numeric(object@meta.data[[m]])
    d <- stats::density(clusid)

    # Map embedding columns to names similar to DimPlot
    colnames_map <- switch(
      reduction,
      tsne = c("tSNE_1","tSNE_2"),
      umap = c("UMAP_1","UMAP_2"),
      pca  = c("PC_1","PC_2")
    )
    coords <- data.frame(
      ident = as.character(idents),
      umap1 = emb[, 1],
      umap2 = emb[, 2],
      clusid = clusid,
      check.names = FALSE
    )
    colnames(coords)[2:3] <- colnames_map

    coords <- dplyr::mutate(coords, sample_clusid = coords$clusid)
    coords <- dplyr::arrange(coords, clusid)

    clusid.df <- data.frame(id = object@meta.data[[group.var]],
                ModuleScore = object@meta.data[[m]],
                stringsAsFactors = FALSE)

    res[[celltype_name]] <- list(
      object = object,
      m = m,
      coords = coords,
      clusid.df = clusid.df,
      density = d
    )
  }

  res
}

#' @export
build_modscore_plots <- function(object, m, coords, clusid.df, d, threshold,
                                 gradient.ft.size = 6, violin.ft.size = 6, step.size = 0.1,
                                 group.var = "orig.ident",
                                 reduction = c("tsne","umap","pca")) {
  reduction <- match.arg(reduction)

  g <- ggplot2::ggplot(coords, ggplot2::aes(x = coords[[2]], y = coords[[3]])) +
    ggplot2::theme_bw() + ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::geom_point(ggplot2::aes(colour = sample_clusid), alpha = 0.5, shape = 20, size = 1) +
    ggplot2::scale_color_gradientn(colours = c("blue4", "lightgrey", "red"),
      values = scales::rescale(c(0, threshold/2, threshold, (threshold + 1)/2, 1), limits = c(0, 1))) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 5, alpha = 1))) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) +
    ggplot2::xlab(if (reduction == "tsne") "tsne-1" else if (reduction == "umap") "umap-1" else "pc-1") +
    ggplot2::ylab(if (reduction == "tsne") "tsne-2" else if (reduction == "umap") "umap-2" else "pc-2")

  g1 <- Seurat::RidgePlot(object, features = m, group.by = group.var) +
    ggplot2::theme(legend.position = "none", title = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = gradient.ft.size)) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "red3") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, step.size))

  g3 <- ggplot2::ggplot(data.frame(x = d$x, y = d$y), ggplot2::aes(x, y)) +
    ggplot2::xlab("ModuleScore") + ggplot2::ylab("Density") + ggplot2::geom_line() +
    ggplot2::geom_segment(ggplot2::aes(xend = d$x, yend = 0, colour = x)) + ggplot2::scale_y_log10() +
    ggplot2::scale_color_gradientn(colours = c("blue4", "lightgrey", "red"),
      values = scales::rescale(c(0, threshold/2, threshold, (threshold + 1)/2, 1), limits = c(0, 1))) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "red3") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, step.size)) + ggplot2::theme(legend.title = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 6))

  arranged <- gridExtra::arrangeGrob(g, g1, g3, ncol = 2, top = grid::textGrob(m, gp = grid::gpar(fontsize = 14, fontface = "bold")))
  list(g = g, g1 = g1, g3 = g3, arranged = arranged)
}
