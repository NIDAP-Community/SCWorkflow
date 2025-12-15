#' @title Update metadata slot of Seurat-class object with custom labels and 
#' provide plot with percentage of cell types
#'  
#' @description Maps custom cluster names to Seurat Object cluster IDs and adds
#'  cluster names to a new metadata column called *Clusternames*. Provides a 
#'  dotplot of percentage of cell types within each cluster. 
#'
#' @param object Seurat-class object with cluster IDs column and cell type 
#'  column present
#' @param cluster.numbers Vector containing cluster numbers that match the 
#'  (numeric) cluster ID's in the cluster.column in Seurat Object metadata
#' @param cluster.names Vector containing custom cluster labels
#' @param cluster.column Column name containing cluster ID in the metadata slot
#'  in the object
#' @param labels.column Column name containing labels (usually cell type) in the
#'  metadata slot in the object
#' @param order.clusters.by Vector containing order of clusters in graph. Can 
#'  contain a subset of cluster numbers to plot that match at least some of
#'  the values in the cluster.column. If NULL, use default order 
#'  (default is NULL)
#' @param order.celltypes.by Vector containing order of cell types in graph. 
#'  Can contain a subset of cell types to plot that match at least some of the
#'  values in the labels.column. If NULL, use default order (default is NULL)
#' @param interactive If TRUE, draw plotly plot (default is FALSE)
#'
#' @importFrom dplyr pull
#' @importFrom reshape2 reshape2::melt
#' @importFrom ggplot2 ggplot geom_point aes theme_classic ylim scale_y_reverse
#' theme ggtitle
#' @importFrom plotly ggplotly
#' @importFrom Seurat AddMetaData
#' @importFrom tibble deframe tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_sort
#'
#' @export
#' @return Returns Seurat-class object with updated meta.data slot containing
#' custom cluster annotation and a plot

nameClusters <- function(object,
                         cluster.identities.table,
                         cluster.numbers,
                         cluster.names,
                         cluster.column,
                         labels.column,
                         order.clusters.by = NULL,
                         order.celltypes.by = NULL,
                         interactive = FALSE)
{
  
  # # Assign cluster numbers with custom cluster names 
  # cluster.numbers <- factor(cluster.numbers,
  #                           str_sort(cluster.numbers, numeric = TRUE))
  # levels(cluster.numbers) <-
  #   str_sort(cluster.numbers, numeric = TRUE)
  # # names(cluster.names) <- as.character(cluster.numbers)
  
  cluster.names= setNames(cluster.identities.table[[cluster.names]],cluster.identities.table[[cluster.numbers]])
  
  
  # get metadata from object and get cluster column
  metadata.df <- object@meta.data
  colval <- metadata.df[[cluster.column]]
  
  
  # # If cluster numbers on input table match the cluster numbers in cluster 
  # # column in metadata, add new custom labels to seurat object
  if (all(unique(metadata.df[[cluster.column]]) %in% unique(cluster.identities.table[[cluster.numbers]]))) {
    object =
    AddMetaData(object,
                metadata =
                  deframe(tibble(metadata.df[[cluster.column]],
                                 cluster.names[as.character(colval)])),
                          col.name = "Clusternames")
  } else{
    stop("Cluster ID's have to match metadata column. Please check entry in
           input table.")
  }
  
  clus.num <-
    as.data.frame.matrix(table(object@meta.data$Clusternames,
                               object@meta.data[[labels.column]]))
  clusnum.df <- reshape2::melt(as.matrix(clus.num))
  sums <- rowSums(clus.num)
  cluster.perc <- (clus.num / sums) * 100
  
  # draw plot 
  clus.df <- reshape2::melt(as.matrix(cluster.perc))
  clus.df$num <- clusnum.df$value
  colnames(clus.df) <-
    c("cluster", "celltype", "percent", "number")
  head(clus.df)
  clus.df$cluster <- factor(clus.df$cluster)
  
  #For plot, optional ordering by clusters
  if (length(order.clusters.by) > 0) {
    add2 <- order.clusters.by[!order.clusters.by %in%
                                levels(clus.df$celltype)]
    if (length(add2) > 0) {
      add2 <- paste(add2, collapse = ", ")
      warning(sprintf("Some factors are not in data: %s", add2))
    }
    clus.dat2 <- unique(clus.df$cluster)
    minus2 <- clus.dat2[!clus.dat2 %in% order.clusters.by]
    if (length(minus2) > 0) {
      minus2 <- paste(minus2, collapse = ", ")
      warning(sprintf("Some factors were not included in the list: %s",
                      minus2))
    }
    clus.df <-
      clus.df %>% filter(cluster %in% order.clusters.by)
    clus.df$cluster <-
      factor(clus.df$cluster, levels = order.clusters.by)
  } else {
    clus.df$cluster <- factor(clus.df$cluster,
                              levels = str_sort(unique(cluster.names), numeric = TRUE))
  }
  
  #For plot, optional if ordering by cell types
  if (length(order.celltypes.by) > 0) {
    #Check for factors not present in dataset
    add <- order.celltypes.by[!order.celltypes.by %in%
                                levels(clus.df$celltype)]
    if (length(add) > 0) {
      add <- paste(add, collapse = ", ")
      warning(sprintf("Some factors are not in data: %s", add))
    }
    clus.dat <- unique(clus.df$celltype)
    minus <- clus.dat[!clus.dat %in% order.celltypes.by]
    if (length(minus) > 0) {
      minus <- paste(minus, collapse = ", ")
      warning(sprintf("Some factors were not included in the list: %s",
                      minus))
    }
    clus.df <- 
      clus.df %>% filter(celltype %in% order.celltypes.by)
    clus.df$celltype <-
      factor(clus.df$celltype, levels = order.celltypes.by)
  }
  
  clus.df <- clus.df %>%
    arrange(cluster, celltype)

  ## if % less then 5 do not plot (convert to NA)
  clus.df[clus.df<5]=NA


  ## print clusters in log:
  clus <- clus.df %>%
    pull(cluster) %>%
    unique() %>%
    as.character()
  cat("clusters:\n")
  cat(clus, sep = "\n")
  cat("\n\n")
  
  ## print cell types in log
  celltypes <- clus.df %>%
    pull(celltype) %>%
    unique() %>%
    as.character()
  cat("celltypes:\n")
  cat(celltypes, sep = "\n")
  
  
  # do plot (suppressMessages for ggplot2 scale replacemnt)
  g <- ggplot(clus.df,
              aes(
                x = celltype,
                y = cluster,
                size = percent,
                color = celltype,
                label = number
              )) +
    theme_classic() +
    geom_point(alpha = 0.5) +
    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )) +
    ggtitle("Percentage of Cell Type within Clusters")
  
  if (interactive == TRUE) {
    g <- ggplotly(g)
  }

  return(list(
    "object" = object,
    "plots" = g
  ))
  
}
