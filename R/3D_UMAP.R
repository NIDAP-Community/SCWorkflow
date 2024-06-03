#' @title Plot 3D-UMAP given a Seurat Object and returns plotly image
#' @description This method provides visualization of 3D-UMAP plot given a
#'   Seurat Object and returns a plotly plot and a dataframe of UMAP
#'   coordinates. It optionally saves the plotly image embedded in an html file.
#'
#' @param object Seurat-class object
#' @param color.variable Metadata column in Seurat Object to use for color
#' @param label.variable Metadata column in Seurat Object to use for label
#' @param dot.size Dot size for plot (default is 4)
#' @param legend If TRUE, show legend (default is TRUE)
#' @param colors Colors used for the color.variable
#' @param filename Filename for saving plot (default is "plot.html")
#' @param npcs Number of principal components used for UMAP calculations
#'   (default is 15)
#' @param save.plot Save plot as widget in html file (default is FALSE)
#'
#' @importFrom plotly as_widget plot_ly
#' @importFrom Seurat RunUMAP
#' @importFrom htmlwidgets saveWidget
#'
#' @export

UMAP3D <- function(object,
                   color.variable,
                   label.variable,
                   dot.size = 4,
                   legend = TRUE,
                   colors = c("darkblue","purple4","green","red","darkcyan",
                              "magenta2","orange","yellow","black"),
                   filename = "plot.html",
                   save.plot = FALSE,
                   npcs = 15) {
  
  #Run UMAP again to get 3d coordinates:
  object <- RunUMAP(
    object,
    assay = "SCT",
    dims = 1:npcs,
    n.components = 3,
    seed.use = 1
  )
  
  umap.coord <-
    as.data.frame(object@reductions$umap@cell.embeddings)
  
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
  umap.df <- data.frame(
    UMAP1 = umap.coord$UMAP_1,
    UMAP2 = umap.coord$UMAP_2,
    UMAP3 = umap.coord$UMAP_3,
    colorvar = as.factor(object@meta.data[[color.variable]]),
    label = object@meta.data[[label.variable]]
  )
  
  fig <- plot_ly(
    umap.df,
    x = ~ UMAP1,
    y = ~ UMAP2,
    z = ~ UMAP3,
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
  
  umap.results <- list("plot"  = fig, "data" = umap.df)
  return(umap.results)
}
