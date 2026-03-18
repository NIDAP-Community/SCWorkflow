# This code comes from NIDAP 'Recluster Filtered Seurat Object'

#' @title Recluster Filtered Seurat Object.
#' @description This template reclusters a filtered Seurat object.
#' @details This method reclusters the filtered SO, preserving the original
#' SCT clustering columns with a prepended prefix, and making new SCT clustering
#' columns based on the reclustering. The image returned is the reclustered
#' project.
#'
#' @param object The input Seurat Object.
#' @param prepend.txt Text to prepend to old columns to make them unique from
#' new. Default is "old".
#' @param old.columns.to.save Old seurat clustering columns (e.g.
#' SCT_snn_res.0.4) to save.
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
#'
#' @return Function returns a reclustered Seurat Object with new clustering
#' columns and renamed original clustering columns, along with a plot of the new
#' dimsensionality reduction.

# Recluster Filtered Seurat Object [scRNA-seq][CCBR]
# (576fe688-c445-48c6-a40a-033da171c149): v11
reclusterFilteredSeuratObject <- function(object,
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
  
  ## Check columns.
  if (length(old.columns.to.save) > 0 &
      !all(old.columns.to.save %in% colnames(object[[]]))) {
    ## ERROR PRESENT ON LINE BELOW. RELATES TO DOTS & UNDERSCORES NIDAPism
    colnames(object[[]]) <- gsub("\\.", "_", colnames(object[[]]))
    if (!all(old.columns.to.save %in% colnames(object[[]]))) {
      stop("Could not find requested metadata columns!")
    }
  }
  
  ## Get columns.
  for (i in seq_along(old.columns.to.save)) {
    old.column.name <- old.columns.to.save[i]
    new.column.name <- paste(prepend.txt, old.column.name, sep = "_")
    object@meta.data[[new.column.name]] <-
      object@meta.data[[old.column.name]]
  }
  
  ## Remove original cluster columns because they are inaccurate
  columns_to_remove <-
    c("seurat_clusters",
      grep("^SCT_snn_res", colnames(object[[]]), value = TRUE))
  object@meta.data %>% select(-one_of(columns_to_remove)) -> object@meta.data
  
  ## Find new clusters.
  object <-
    FindVariableFeatures(
      object = object,
      nfeatures = 2000,
      mean.cutoff = c(0.1, 8),
      dispersion.cutoff = c(1, 100000),
      selection.method = "vst",
      verbose = FALSE
    )
  
  # Method for automatically picking NPCS
  if (number.of.pcs == 0) {
    object <-
      RunPCA(
        object = object,
        npcs = 30,
        verbose = FALSE,
        seed.use = 42
      ) # initial run
    pca_stdev <- Seurat::Stdev(object, reduction = "pca")
    sumpcsd = sum(pca_stdev)
    pcvar = (pca_stdev / sumpcsd) * 100
    cumu <- cumsum(pcvar)
    co1 <- which(cumu > 80 & pcvar < 5)[1]
    co2 <-
      sort(which((pcvar[1:length(pcvar) - 1] - pcvar[2:length(pcvar)]) > 0.1), 
           decreasing = T)[1] + 1
    number.of.pcs = min(co1, co2)
    print(number.of.pcs)
  }
  
  ## Dimensionality reduction
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
      dims = 1:number.of.pcs,
      seed.use = 1
    )
  object <- FindNeighbors(object, dims = 1:number.of.pcs)
  
  ## Find Clusters.
  resolutions <-
    seq(
      cluster.resolution.low.range,
      cluster.resolution.high.range,
      cluster.resolution.range.bins
    )
  
  for (r in resolutions) {
    object <- FindClusters(object, resolution = r, algorithm = 1)
  }
  print("Clustering successful!")
  
  ## Fix orig_ident back to orig.ident.
  colnames(object@meta.data)[colnames(object@meta.data) == "orig_ident"] <-
    "orig.ident"
  
  plot.list <-
    lapply(paste0("SCT_snn_res.", resolutions), function(z)
      DimPlot(object, reduction = reduction.type, group.by = z) +
        labs(title = z) +
        theme(plot.title = element_text(hjust = 0.5)))
  
  ## Not sure if we need this block with 'g'.
  ncol = ceiling(length(resolutions) ^ 0.5)
  nrow = ceiling(length(resolutions) / ncol)
  g <- plot_grid(plotlist = plot.list,
                 nrow = nrow,
                 ncol = ncol)
  # print(g)
  
  # return(list("object"=object, "plot"=plot.list))
  return(
    list(
      "object" = object, 
      "plots" = g
        )
      )
}
