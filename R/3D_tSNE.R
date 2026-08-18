#' @title Plot 3D-TSNE given a Seurat Object and returns plotly image
#' @description This method provides visualization of a 3D t-SNE or UMAP plot
#'   given a Seurat Object and returns a plotly plot and a dataframe of embedding
#'   coordinates. It optionally saves the plotly image embedded in an html file.
#'
#' @param object Seurat-class object
#' @param color.variable Metadata column in Seurat Object to use for color
#' @param label.variable Metadata column in Seurat Object to use for label
#' @param plot.type Dimensionality reduction to calculate and plot. One of
#'   "tSNE" or "UMAP" (default is "tSNE")
#' @param dot.size Dot size for plot (default is 4)
#' @param legend If TRUE, show legend (default is TRUE)
#' @param colors Colors used for the color.variable
#' @param filename Filename for saving plot (default is "plot.html")
#' @param npcs Number of principal components used for t-SNE or UMAP calculations
#'   (default is 15)
#' @param save.plot Save plot as widget in html file (default is FALSE)
#'
#' @importFrom plotly as_widget plot_ly
#' @importFrom Seurat RunTSNE RunUMAP
#' @importFrom htmlwidgets saveWidget
#'
#' @export
#'
#' @return A list with a plotly 3D TSNE plot (`figure`) and TSNE coordinates
#' (`tsne.df`).
#'
#' @examples
#' \dontrun{
#' out <- tSNE3D(
#'   object = seurat_obj,
#'   color.variable = "cell_type",
#'   label.variable = "orig.ident",
#'   plot.type = "UMAP",
#'   npcs = 15,
#'   save.plot = FALSE
#' )
#' }

tSNE3D <- function(object,
                   color.variable,
                   label.variable,
                   plot.type = "tSNE",
                   dot.size = 4,
                   legend = TRUE,
                   colors = c("darkblue","purple4","green","red","darkcyan",
                              "magenta2","orange","yellow","black"),
                   filename = "plot.html",
                   save.plot = FALSE,
                   npcs = 15) {
  plot.type <- trimws(as.character(plot.type))
  plot.type.normalized <- tolower(plot.type)
  if (!plot.type.normalized %in% c("tsne", "umap")) {
    stop("plot.type must be one of: tSNE, UMAP", call. = FALSE)
  }

  if (identical(plot.type.normalized, "umap")) {
    object <- RunUMAP(
      object,
      assay = "SCT",
      dims = 1:npcs,
      n.components = 3,
      seed.use = 1
    )
    embedding.coord <- as.data.frame(object@reductions$umap@cell.embeddings)
    output.columns <- c("UMAP_1", "UMAP_2", "UMAP_3")
  } else {
    object <- RunTSNE(
      object,
      assay = "SCT",
      dims = 1:npcs,
      dim.embed = 3,
      seed.use = 1
    )
    embedding.coord <- as.data.frame(object@reductions$tsne@cell.embeddings)
    output.columns <- c("tSNE_1", "tSNE_2", "tSNE_3")
  }
  
  if (is.null(object@meta.data[[label.variable]])) {
    stop(
      sprintf(
        "The metadata variable selected for labeling %s is not
                 available in the seurat object",
        label.variable
      )
    )
  }
  if (is.null(object@meta.data[[color.variable]])) {
    stop(
      sprintf(
        "The metadata variable selected for color %s is not available
                 in the seurat object",
        color.variable
      )
    )
  }
  
  #Set up dataframe for plotly:
  embedding.df <- data.frame(
    dim1 = embedding.coord[[1]],
    dim2 = embedding.coord[[2]],
    dim3 = embedding.coord[[3]],
    colorvar = as.factor(object@meta.data[[color.variable]]),
    label = object@meta.data[[label.variable]]
  )
  colnames(embedding.df)[1:3] <- output.columns
  
  fig <- plot_ly(
    embedding.df,
    x = stats::as.formula(paste0("~", output.columns[[1]])),
    y = stats::as.formula(paste0("~", output.columns[[2]])),
    z = stats::as.formula(paste0("~", output.columns[[3]])),
    color = ~ colorvar,
    colors = colors,
    type = "scatter3d",
    mode = "markers",
    hoverinfo = 'text',
    text = ~ label,
    size = dot.size
  )
  
  if (legend == FALSE) {
    fig <- hide_legend(fig)
  }
  
  #Save plot into html as embedded plotly image
  if (save.plot == TRUE) {
    htmlwidgets::saveWidget(as_widget(fig), filename, selfcontained = TRUE)
  }
  
  tsne.results <- list( "data" = embedding.df,
                        "plots"  = fig)
  return(tsne.results)
}
