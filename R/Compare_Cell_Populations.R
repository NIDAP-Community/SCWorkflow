#' @title Compare Cell Populations
#' @description Compare cell population distributions across different groups
#' using bar plots and box plots. Creates visualizations showing cell type
#' frequencies or counts across user-defined groupings.
#' 
#' @details This function generates comparative visualizations of cell 
#' populations from a Seurat object. It can display data as either frequency 
#' percentages or absolute counts, and creates both stacked bar plots 
#' (with alluvial flow connections) and grouped box plots for comparison 
#' across samples and conditions.
#'
#' @param object A Seurat object containing the single-cell data
#' @param annotation.column Character string specifying the metadata column 
#' containing cell type annotations to summarize in the bar plot
#' @param group.column Character string specifying the metadata column 
#' defining groups to compare (e.g., treatment conditions)
#' @param sample.column Character string specifying the metadata column 
#' containing sample identifiers. Default is "orig.ident"
#' @param counts.type Character string specifying plot data type: 
#' "Frequency" (percentages) or "Counts" (absolute numbers). Default is "Frequency"
#' @param group.order Character vector specifying the order of groups in plots.
#' If NULL, uses natural order from data. Default is NULL
#' @param seurat.object.filename Character string for the Seurat object 
#' filename. Default is "seurat_object.rds"
#' @param wrap.ncols Integer specifying number of columns for facet wrapping 
#' in box plots. Default is 5
#'
#' @import Seurat
#' @import ggplot2
#' @import ggpubr
#' @import RColorBrewer
#' @import tibble
#' @import reshape2
#' @import data.table
#' @import dplyr
#' @import magrittr
#' @import cowplot
#' @import gridExtra
#' @import grid
#' @import scales
#'
#' @importFrom ggalluvial geom_flow
#' @importFrom stats setNames
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Plots} - A list with two ggplot objects:
#'     \itemize{
#'       \item \code{Barplot} - Stacked bar plot with alluvial flows
#'       \item \code{Boxplot} - Faceted box plots by cell type (only if counts.type="Frequency")
#'     }
#'   \item \code{Table} - A data.frame with cell counts and percentages
#' }
#'
#' @examples
#' \dontrun{
#' # Compare cell populations by treatment group
#' results <- compareCellPopulations(
#'   object = seurat_obj,
#'   annotation.column = "cell_type",
#'   group.column = "treatment",
#'   sample.column = "sample_id",
#'   counts.type = "Frequency"
#' )
#' 
#' # Display plots
#' plot(results$Plots$Barplot)
#' plot(results$Plots$Boxplot)
#' 
#' # View summary table
#' head(results$Table)
#' }

compareCellPopulations <- function(
  object,
  annotation.column,
  group.column,
  sample.column = "orig.ident",
  counts.type = "Frequency",
  group.order = NULL,
  wrap.ncols = 5
) {
  
  ## -------------------------------- ##
  ## Input Validation                 ##
  ## -------------------------------- ##
  
  # Validate object
  if (!inherits(object, "Seurat")) {
    stop("Error: 'object' must be a Seurat object")
  }
  
  # Validate counts.type
  if (!counts.type %in% c("Frequency", "Counts")) {
    stop("Error: 'counts.type' must be either 'Frequency' or 'Counts'")
  }
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  createAnnoTable <- function(SO, AnnoCol, GroupCol) {
    ## Extract annotation data for each group using a 2D contingency table
    cntMat <- table(SO@meta.data[[AnnoCol]], SO@meta.data[[GroupCol]])
    
    # Convert to data frame while preserving row/column names
    cntTble <- as.data.frame.matrix(cntMat)
    cntTble <- data.frame(
      lapply(cntTble, function(x) as.numeric(as.character(x))),
      check.names = FALSE,
      row.names = rownames(cntTble)
    )
    
    freqTble <- apply(cntTble, 2, FUN = function(x) {
      return(x / sum(x))
    })
    freqTble <- (freqTble * 100)
    
    outTbl <- merge(cntTble, as.data.frame(freqTble),
                    by = 'row.names',
                    suffixes = c('_CellCounts', '_Percent'))
    outTbl <- dplyr::rename(outTbl, 'Clusters' = "Row.names")
    
    return(list(
      'CellFreq' = freqTble,
      'CellCounts' = cntTble,
      'OutTable' = outTbl
    ))
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##

  # Replace dots with underscores in column names
  colnames(object@meta.data) <- gsub("\\.", "_", colnames(object@meta.data))
  
  # Update column names if they were modified
  annotation.column <- gsub("\\.", "_", annotation.column)
  group.column <- gsub("\\.", "_", group.column)
  sample.column <- gsub("\\.", "_", sample.column)
  
  
  # Validate metadata columns exist
  required.cols <- c(annotation.column, group.column, sample.column)
  missing.cols <- setdiff(required.cols, colnames(object@meta.data))
  if (length(missing.cols) > 0) {
    stop("Error: The following columns are missing from metadata: ",
         paste(missing.cols, collapse = ", "))
  }
  
  
  
  
  # Set up ordering
  ordr <- object@meta.data[[annotation.column]] %>% 
    unique() %>% 
    sort()
  
  if (is.null(group.order)) {
    group.order <- unique(object@meta.data[[group.column]])
  }
  
  # Set up colors
  numColors <- max(
    length(unique(object@meta.data[[annotation.column]])),
    20
  )
  colpaired <- colorRampPalette(brewer.pal(12, "Paired"))
  cols <- c(
    "#e6194B", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#42d4f4",
    "#f032e6", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324",
    "#800000", "#aaffc3", "#808000", "#000075",
    colpaired(numColors)
  )
  names(cols) <- ordr
  
  object@meta.data[[annotation.column]] <- factor(
    object@meta.data[[annotation.column]],
    levels = ordr
  )
  
  # Create tables
  ColTables <- createAnnoTable(object, annotation.column, group.column)
  BoxTables <- createAnnoTable(object, annotation.column, sample.column)
  
  metaGroups <- object@meta.data[, c(group.column, sample.column)]
  rownames(metaGroups) <- NULL
  metaGroups <- metaGroups %>% unique()
  
  ## Create plots based on counts type
  if (counts.type == 'Frequency') {
    ptbl <- melt(ColTables$CellFreq)
    ptblBox <- melt(as.matrix(BoxTables$CellFreq))
    ptblBox <- merge(ptblBox, metaGroups, 
                     by.x = 'Var2', by.y = sample.column, all.x = TRUE)
    
    labelCol <- 'PerValue'
    ylab <- 'Frequency of each cell type (100%)'
  } else if (counts.type == "Counts") {
    ptbl <- melt(as.matrix(ColTables$CellCounts))
    ptblBox <- melt(as.matrix(BoxTables$CellCounts))
    ptblBox <- merge(ptblBox, metaGroups, 
                     by.x = 'Var2', by.y = sample.column, all.x = TRUE)
    
    labelCol <- 'value'
    ylab <- 'Cell Counts'
  }
  
  # Format bar plot data
  ptbl$Var1 <- factor(ptbl$Var1, levels = ordr)
  ptbl$value <- round(ptbl$value, 1)
  ptbl$PerValue <- paste0(ptbl$value, '%')
  ptbl$PerValue <- gsub('^%$', "_", ptbl$PerValue)
  ptbl[ptbl$value < 1, 'PerValue'] <- ""
  ptbl$Var2 <- factor(ptbl$Var2, levels = group.order)
  
  # Create bar plot with alluvial flows
  p2 <- ptbl %>%
    ggplot(aes_string(y = 'value', x = 'Var2', fill = 'Var1', label = labelCol)) +
    geom_flow(aes(alluvium = Var1), alpha = .2,
              lty = 2, color = "black",
              curve_type = "linear",
              width = .5) +
    geom_col(aes(fill = Var1), width = .5, color = "black") +
    geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    theme_classic() +
    ylab(ylab) +
    xlab("") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_manual(annotation.column, values = cols)
  
  # Create box plot
  ptblBox$value <- round(ptblBox$value, 1)
  ptblBox$PerValue <- paste0(ptblBox$value, '%')
  ptblBox$PerValue <- gsub('^%$', "_", ptblBox$PerValue)
  ptblBox[ptblBox$value < 1, 'PerValue'] <- ""
  ptblBox[[group.column]] <- factor(ptblBox[[group.column]], levels = group.order)
  
  p2_Box <- ptblBox %>%
    ggboxplot(y = 'value', x = group.column, add = "jitter", color = "Var1") +
    facet_wrap(~Var1, ncol = wrap.ncols, scales = 'fixed') +
    ylab(ylab) +
    xlab("") +
    theme(legend.title = element_blank())
  
  # Return results
  result <- list(
    'plots' = list('Barplot' = p2, 'Boxplot' = p2_Box),
    'data' = ColTables$OutTable
  )
  
  return(result)
}

# Add global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
  "Var1", "Var2", "value", "PerValue", "alluvium",
  ".", "CellFreq", "CellCounts", "OutTable"
))
