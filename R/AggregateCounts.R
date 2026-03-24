#' @title Aggregate Counts (Pseudobulk)
#' @description Compute pseudobulk expression by averaging expression across groups
#'              defined by one or more metadata columns, and return a tidy table.
#' @details Uses Seurat's `AverageExpression()` on the `SCT` assay to compute
#'          group-wise average expression for each feature. Also produces a
#'          bar plot (via `ggplot2`/`plotly`) showing the number of cells per
#'          pseudobulk group and warns if any group contains only one cell.
#'
#' @param object Seurat-class object.
#' @param var.group Character vector of metadata column names used to define
#'                  pseudobulk groups. When multiple columns are supplied, an
#'                  interaction of these columns defines the groups.
#' @param slot Character name of the assay data layer passed to
#'             `AverageExpression()` (e.g., "data", "counts", or "scale.data").
#'
#' @import Seurat
#' @import tidyverse
#' @import ggplot2
#' @import plotly
#' @importFrom dplyr select
#'
#' @export
#'
#' @return A data.frame of pseudobulk expression with columns `Gene` followed by
#'         one column per pseudobulk group. Column names are sanitized to
#'         contain only alphanumeric/underscore characters.
#'
#' @examples
#' \dontrun{
#' out <- aggregateCounts(
#'   object = seurat_obj,
#'   var.group = c("orig.ident", "condition"),
#'   slot = "data"
#' )
#' }

aggregateCounts <- function(object,
                            var.group,
                            slot){
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  pseudobulk <- AverageExpression(object, 
                                  return.seurat = FALSE,
                                  assay = "SCT",
                                  group.by = var.group,
                                  layer = slot)[[1]] %>% 
                  as.data.frame.matrix() 
    
  pseudobulk$Gene <- rownames(pseudobulk)
  pseudobulk <- pseudobulk %>% select("Gene", everything())
  rownames(pseudobulk) <- NULL

  # Further processing of column names
  colnames(pseudobulk) <- gsub("\\W","_",colnames(pseudobulk))
  
  # Return Table/Figure that gives statistics on 
  # Number of Cells in each group/new sample
  # Distribution of cell Counts in each group/new sample
  meta <- object@meta.data[,var.group]
  
  # check that columns are all factors / categorical
  char_or_factor_cols <- sapply(meta, function(x) is.character(x) || is.factor(x) || is.logical(x))
  
  # do plots and tables if all columns are factors or characters
  if(all(char_or_factor_cols)){
    meta$interaction <- gsub("\\W","_",interaction(meta))
    
    df <- as.data.frame(table(pseudobulk_group = meta$interaction)) %>% filter(Freq != 0)
    # sort the table by the number of cells in each group
    df <- df[order(df$Freq, decreasing = F),]
    
    if(any(df$Freq == 1)){
      single_counts <- df$pseudobulk_group[df$Freq == 1]
      # sprintf, make custom warning message with %s as placeholder for single count groups
      warning(sprintf(
        "Some groups have only 1 cell. It is recommended to have at least 2 cells in each group.\nAffected groups: %s",
        paste(single_counts, collapse = ", ")
      ))
    }
    
    p <- ggplotly(ggplot(df, aes(x = pseudobulk_group, y = Freq)) +
                    geom_bar(stat = "identity", position = "stack") +
                    labs(y = "Counts", x = "Pseudobulk Groups", title = "Number of Cells in each Pseudobulk Group") +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    
    
  } else {
    stop("All columns in var.group must be factors or characters")
  }
  
  return(list(data=pseudobulk,
              plot=p))
}