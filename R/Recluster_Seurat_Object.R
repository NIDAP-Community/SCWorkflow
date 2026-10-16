
#' @title Recluster Seurat Object.
#' @description This template reclusters a Seurat object.
#' @details This method reclusters the input SO, preserving the original SCT 
#' clustering columns with a prepended prefix, and making new SCT clustering
#' columns based on the reclustering. The image returned is the reclustered
#' project.
#'
#' @param object The input Seurat Object.
#' @param prepend.txt Text to prepend to old columns to make them unique from
#' new. Default is "old".
#' @param old.columns.to.save Old seurat clustering columns (e.g.
#' SCT_snn_res.0.4) to save. Select only columns that will be replaced during
#' reclustering.
#' @param number.of.pcs Select the number of principal components for your
#' analysis. Set to 0 to automatically decide. Default is 50.
#' @param cluster.resolution.low.range Select minimum resolution for clustering
#' plots. The lower you set this, the FEWER clusters will be generated. Default
#' is 0.2.
#' @param cluster.resolution.high.range Select maximum resolution for clustering
#' plots. The higher you set this, the MORE clusters will be generated. Default
#' is 1.2.
#' @param cluster.resolution.range.bins Select the bins for your cluster plots.
#' For example, if you input 0.2 as your bin, and have low/high resolution
#' ranges of 0.2 and 0.6, then the template will produce cluster plots at
#' resolutions of 0.2, 0.4 and 0.6. Default is 0.2.
#' @param reduction.type  Select the kind of clustering visualization you would
#' like to use to visualize the cell type results ("umap", "tsne", "pca").
#' Default is "tsne".
#'
#' @import Seurat
#' @import cowplot
#' @import tidyverse
#'
#' @export
#' @return Function returns a reclustered Seurat Object with new clustering
#' columns and renamed original clustering columns, along with a plot of the new
#' dimsensionality reduction.

## Origin: This code adapted from v11 of our Potomac released NIDAP template:
## Recluster Filtered Seurat Object [scRNA-seq][CCBR]
## (576fe688-c445-48c6-a40a-033da171c149): v11

reclusterSeuratObject <- function(object,
                                  prepend.txt = "old",
                                  old.columns.to.save,
                                  number.of.pcs = 50,
                                  cluster.resolution.low.range = 0.2,
                                  cluster.resolution.high.range = 1.2,
                                  cluster.resolution.range.bins = 0.2,
                                  reduction.type = "tsne"
                                  ) {

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##

  if (length(number.of.pcs) != 1 ||
      !is.numeric(number.of.pcs) ||
      !is.finite(number.of.pcs) ||
      number.of.pcs < 0 ||
      number.of.pcs != floor(number.of.pcs)) {
    stop("number.of.pcs must be a non-negative whole number.")
  }

  resolution.parameters <- list(
    cluster.resolution.low.range,
    cluster.resolution.high.range,
    cluster.resolution.range.bins
  )
  if (!all(vapply(
    resolution.parameters,
    function(value) is.numeric(value) && length(value) == 1 && is.finite(value),
    logical(1)
  ))) {
    stop("Cluster resolution parameters must be finite numeric values.")
  }
  if (cluster.resolution.low.range > cluster.resolution.high.range) {
    stop(
      "cluster.resolution.low.range must be less than or equal to ",
      "cluster.resolution.high.range."
    )
  }
  if (cluster.resolution.range.bins <= 0) {
    stop("cluster.resolution.range.bins must be greater than zero.")
  }

  metadata.column.names <- colnames(object[[]])
  clustering.graph.name <- paste0(DefaultAssay(object), "_snn")
  assay.graph.names <- as.vector(outer(
    names(object@assays),
    c("_nn", "_snn"),
    paste0
  ))
  clustering.graph.names <- unique(c(
    names(object@graphs),
    assay.graph.names
  ))
  clustering.column.prefix <- paste0(clustering.graph.name, "_res.")
  normalized.clustering.prefixes <- paste0(
    gsub("[^[:alnum:]]", "", clustering.graph.names),
    "res"
  )
  normalized.metadata.column.names <- gsub(
    "[^[:alnum:]]",
    "",
    metadata.column.names
  )
  clustering.columns.to.replace <- unique(c(
    intersect(
      c("seurat_clusters", "seurat.clusters"),
      metadata.column.names
    ),
    metadata.column.names[vapply(
      normalized.metadata.column.names,
      function(name) any(startsWith(name, normalized.clustering.prefixes)),
      logical(1)
    )]
  ))
  possible.metadata.columns <- clustering.columns.to.replace
  possible.metadata.message <- paste0(
    "Possible clustering columns to save before replacement:\n",
    paste0("  - ", possible.metadata.columns, collapse = "\n")
  )
  if (length(old.columns.to.save) == 0) {
    message("No old columns were saved.")
  } else {
    missing.metadata.columns <- setdiff(
      old.columns.to.save,
      metadata.column.names
    )
    non.clustering.columns <- setdiff(
      intersect(old.columns.to.save, metadata.column.names),
      possible.metadata.columns
    )
    rejected.columns <- c(missing.metadata.columns, non.clustering.columns)
    if (length(rejected.columns) > 0) {
      rejected.messages <- c()
      if (length(missing.metadata.columns) > 0) {
        rejected.messages <- c(
          rejected.messages,
          paste0(
            "Old column(s) to save not found in metadata: ",
            paste(shQuote(missing.metadata.columns), collapse = ", "),
            "."
          )
        )
      }
      if (length(non.clustering.columns) > 0) {
        rejected.messages <- c(
          rejected.messages,
          paste0(
            "Old column(s) to save are not clustering columns replaced during ",
            "reclustering: ",
            paste(shQuote(non.clustering.columns), collapse = ", "),
            ".\n     These metadata columns are retained automatically."
          )
        )
      }
      rejected.columns.message <- paste(rejected.messages, collapse = "\n")
      old.columns.to.save <- setdiff(old.columns.to.save, rejected.columns)
      if (length(old.columns.to.save) == 0) {
        stop(rejected.columns.message, "\n", possible.metadata.message)
      }
      message(rejected.columns.message)
    }
  }
  
  ## Rename selected clustering columns before creating replacements.
  if (length(old.columns.to.save) > 0) {
    new.column.names <- paste(prepend.txt, old.columns.to.save, sep = "_")
    message(
      "Renaming selected clustering columns before reclustering:\n",
      paste0("  - ", old.columns.to.save, " -> ", new.column.names,
             collapse = "\n")
    )
    column.indices <- match(old.columns.to.save, colnames(object@meta.data))
    colnames(object@meta.data)[column.indices] <- new.column.names
  }
  
  ## Remove original cluster columns because they are about to be replaced.
  columns_to_remove <- setdiff(
    clustering.columns.to.replace,
    old.columns.to.save
  )
  if (length(columns_to_remove) > 0) {
    message(
      "Removing original clustering columns before reclustering:\n",
      paste0("  - ", columns_to_remove, collapse = "\n")
    )
    object@meta.data %>% select(-one_of(columns_to_remove)) -> object@meta.data
  }
  
  ## Find the new reclustered clusters.
  object <-
    FindVariableFeatures(
      object = object,
      nfeatures = 2000,
      mean.cutoff = c(0.1, 8),
      dispersion.cutoff = c(1, 100000),
      selection.method = "vst",
      verbose = FALSE
    )
  
  # Method for automatically picking NPCs.
  if (number.of.pcs == 0) {
    object <-
      RunPCA(
        object = object,
        npcs = 30,
        verbose = FALSE,
        seed.use = 42
      ) # initial run
    sumpcsd = sum(object@reductions$pca@stdev)
    pcvar = (object@reductions$pca@stdev / sumpcsd) * 100
    cumu <- cumsum(pcvar)
    co1 <- which(cumu > 80 & pcvar < 5)[1]
    co2 <-
      sort(which((pcvar[1:length(pcvar) - 1] - pcvar[2:length(pcvar)]) > 0.1), 
           decreasing = T)[1] + 1
    number.of.pcs = min(co1, co2)
    print(number.of.pcs)
  }
  
  ## Dimensionality reduction.
  object <-
    RunPCA(
      object = object,
      npcs = number.of.pcs,
      verbose = FALSE,
      seed.use = 42
    )
  object <-
    RunUMAP(
      object = object,
      reduction = "pca",
      dims = 1:number.of.pcs,
      seed.use = 42
    )
  object <-
    RunTSNE(
      object = object,
      reduction = "pca",
      dim.embed = 2,
      dims = 1:number.of.pcs,
      ## TJM: Why is this seed different from the one for PCA or UMAP? ====
      seed.use = 1
    )
  object <- FindNeighbors(
    object,
    dims = 1:number.of.pcs
  )
  
  ## Find Clusters.
  resolutions <-
    seq(
      cluster.resolution.low.range,
      cluster.resolution.high.range,
      cluster.resolution.range.bins
    )
  
  for (r in resolutions) {
    object <- FindClusters(
      object,
      resolution = r,
      algorithm = 1
    )
  }
  print("Clustering successful!")
  
  ## Fix orig_ident back to orig.ident.
  colnames(object@meta.data)[colnames(object@meta.data) == "orig_ident"] <-
    "orig.ident"
  
  plot.list <-
    lapply(paste0(clustering.column.prefix, resolutions), function(z)
      DimPlot(object, reduction = reduction.type, group.by = z) +
        labs(title = z) +
        theme(plot.title = element_text(hjust = 0.5)))
  
  ## TJM: Not sure if we need this block with 'g'. ====
  ncol = ceiling(length(resolutions) ^ 0.5)
  nrow = ceiling(length(resolutions) / ncol)
  g <- plot_grid(plotlist = plot.list,
                 nrow = nrow,
                 ncol = ncol)
  ## TJM: Dev output. Remove before release. ====
  # print(g)
  
  # return(list("object"=object, "plot"=plot.list))
  return(
    list(
      "object" = object, 
      "plots" = g
      )
      )
}
