# This code comes from NIDAP 'Dotplot of Gene Expression by Metadata'

#' @title Dotplot of Gene Expression by Metadata
#' @description This function uses the Dotplot function from Seurat and plots
#'        average gene expression values and percent expressed for a set of
#'        genes.
#' @details This method provides a dotplot showing the percent frequency of
#'        gene-positive cells as size of dot and degree of expression as color
#'        of dot.
#'
#' @param object Seurat Object
#' @param metadata Metadata column in Seurat Object to plot
#' @param cells Vector of metadata category factors to plot and should be found
#'  in metadata column. Order of plotting will follow exact order as entered. 
#' @param markers Vector of genes to plot.  Order of plotting will follow exact
#'  order as entered
#' @param plot.reverse If TRUE, set metadata categories to x-axis and genes to 
#'  y-axis (default is FALSE)
#' @param cell.reverse.sort If TRUE, Reverse plot order of metadata category 
#'  factors (default is FALSE)
#' @param dot.color Dot color (default is "dark blue")
#' @importFrom tidyr pivot_wider
#' @importFrom Seurat Idents DotPlot
#' @export 
#' 
#' @return Dotplot with markers and cell types. 

dotPlotMet <- function(object,
                       metadata,
                       cells,
                       markers,
                       plot.reverse = FALSE,
                       cell.reverse.sort = FALSE,
                       dot.color = "darkblue") {

### ignore extra labels:
metadata.df <- object@meta.data
ExtraValue <- sum(!cells %in% unique(metadata.df[[metadata]]))
if (ExtraValue > 0) {
         missinglab2 <- cells[!cells %in% unique(metadata.df[[metadata]])]
         warning(sprintf("There are %s additional elements in your input categories\n that are missing from your metadata table: ", ExtraValue))
         missinglab2 <- cat(paste(as.character(missinglab2), collapse = "\n"))
#
    cells <- cells[cells %in% unique(metadata.df[[metadata]])]
#
    }
### End of the check  
  
  #Set up metadata as new identity:
  metadata.df <- object@meta.data
  Idents(object) <- metadata.df[[metadata]]
  
  #### Error messages #### 
  
  #Calculate difference in input category elements from metadata column elements
  a <- sum(cells %in% unique(Idents(object)))
  b <- sum(!unique(Idents(object)) %in% cells)
  c <- sum(!cells %in% unique(Idents(object)))
  
  if (a < 2) {
    stop(
      "At least 2 metadata categories you wish to plot should
         match the metadata column you are plotting."
    )
  }
  
  if (b > 0) {
    missinglab <-
      unique(Idents(object))[!unique(Idents(object)) %in% cells]
    warning(
      sprintf(
        "There are %s additional element(s) in the metadata table category that 
        were missing from your input categories: %s", b,
        paste(as.character(missinglab), collapse = "\n")
      )
    )
  }
  
  if (c > 0) {
    missinglab2 <- cells[!cells %in% unique(Idents(object))]
    warning(
      sprintf(
        "There are %s additional elements in your input categories
          that are missing from your metadata table: ", c)
      )
    missinglab2 <-
      cat(paste(as.character(missinglab2), collapse = "\n"))
  }
  
  #Subset object by identity
  object <- subset(object, idents = cells)
  
  #Clean up input genes and custom names
  markers <- gsub("[[:space:]]", "", markers)
  l1 <- length(markers)
  dups <- markers[duplicated(markers)]
  markers <- markers[!duplicated(markers)]
  
  l2 <- length(markers)
  print(sprintf("There are %s total unique genes in the genelist", l2))
  if (l1 > l2) {
    warning(sprintf("\n\nThe following duplicate genes were removed: %s",
                   dups))
  }
  missing.genes <- markers[!markers %in% rownames(object)]
  
  l3 <- length(missing.genes)
  missing.genes <- paste(shQuote(missing.genes), collapse = ", ")
  if (l3 == l2) {
    stop("No genes listed are found in dataset.")
  }
  if (l3 > 0) {
    warning(
      sprintf(
        "There are %s gene(s) absent from dataset: %s. Possible reasons are that
          gene is not official gene symbol or gene is not highly expressed and
          has been filtered.\n ", l3, missing.genes) 
      )
  }
  
  #### Main Code: #####
  
  markers <- markers[markers %in% rownames(object)]
  
  cells <- cells[cells != ""]
  
  #Run Seurat Dotplot function
  dp <- DotPlot(object,
                assay = "SCT",
                features = markers,
                dot.scale = 4,
                dot.min = .1,
                cols = c("lightgrey", dot.color)
                )
  cells <- cells[cells %in% dp$data$id]
  if (length(cells) < 1) {
    stop("None of the cell types in list match cell types in metadata")
  }
  cells.missing <- cells[!cells %in% dp$data$id]
  if (length(cells.missing) > 0) {
    cells.missing <- paste(shQuote(cells.missing), collapse = ", ")
    warning(sprintf(
      "Some categories are missing from your dataset: %s",
      cells.missing
    ))
  }
  
  #Add to ggplot object data created by DotPlot and draw plot
  dp$data$id <- factor(dp$data$id, levels = rev(cells))
  dp$data$features.plot <- factor(dp$data$features.plot, levels = markers)
  
  plot <- ggplot(data = dp$data,
                 mapping = aes_string(x = "features.plot", y = "id")) +
            geom_point(mapping = aes_string(size = "pct.exp",
                                    color = "avg.exp.scaled")) +
            scale_color_gradient(low = "lightgrey", high = dot.color) +
            theme_cowplot() +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_text(angle = 90, 
                                             vjust = 1,
                                             hjust = 1)) +
            labs(y = metadata)
  
  if (plot.reverse == TRUE) {
    plot <- plot + coord_flip()
  }
  if (cell.reverse.sort == TRUE) {
    plot <- plot + scale_y_discrete(limits = rev(levels(dp$data$id)))
  }
  
  #Provide Tabular format of Dotplot data 
  dp.pct.tab <- dp$data %>%
    plotly::select(features.plot, pct.exp, id) %>%
    tidyr::pivot_wider(names_from = features.plot,
                       values_from = pct.exp)
  dp.exp.tab <- dp$data %>%
    plotly::select(features.plot, avg.exp.scaled, id) %>%
    tidyr::pivot_wider(names_from = features.plot,
                       values_from = avg.exp.scaled)
  
  result.list <-
    list("data"=list(
            "pct" = dp.pct.tab,
            "exp" = dp.exp.tab),
          "plots" = plot
    )
  
  return(result.list)
  
}
