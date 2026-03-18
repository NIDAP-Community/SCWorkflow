#' @title Add External Cell Annotations
#' @description Merge per-cell external metadata into a Seurat object's metadata
#' using a barcode column.
#'
#' @param object A Seurat object.
#' @param external_metadata A data.frame/data.table containing at least one
#' barcode column and one annotation column.
#' @param barcode_column Column name in external_metadata containing cell
#' barcodes. Default is "Barcode".
#' @param external_cols_to_add Character vector of external columns to import.
#' If NULL, all columns except barcode_column are imported.
#' @param col_to_viz Optional metadata column to visualize on UMAP/tSNE when
#' available.
#'
#' @import Seurat
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggtitle
#' @export
#'
#' @return A list with entries:
#' - object: Seurat object with appended metadata columns
#' - plots: list containing UMAP/tSNE plots when requested and available

AddExternalAnnotation <- function(object,
                                  external_metadata,
                                  barcode_column = "Barcode",
                                  external_cols_to_add = NULL,
                                  col_to_viz = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }

  ext <- as.data.frame(external_metadata)
  if (!(barcode_column %in% colnames(ext))) {
    stop(paste("External metadata does not have the specified barcode column:", barcode_column))
  }

  # Normalize external column selection so both c("A","B") and "A,B" work.
  if (is.null(external_cols_to_add)) {
    cols_to_add <- setdiff(colnames(ext), barcode_column)
  } else if (length(external_cols_to_add) == 1 && grepl(",", external_cols_to_add)) {
    cols_to_add <- trimws(unlist(strsplit(external_cols_to_add, ",")))
  } else {
    cols_to_add <- external_cols_to_add
  }

  missing_cols <- setdiff(cols_to_add, colnames(ext))
  if (length(missing_cols) > 0) {
    stop(paste("These columns were not found in external_metadata:", paste(missing_cols, collapse = ", ")))
  }

  if (length(cols_to_add) == 0) {
    stop("No columns selected to append from external_metadata.")
  }

  meta <- object[[]]
  meta$Barcode <- rownames(meta)

  # Keep first record for duplicated barcodes to avoid row inflation.
  ext <- ext[!duplicated(ext[[barcode_column]]), c(barcode_column, cols_to_add), drop = FALSE]

  matched <- match(meta$Barcode, ext[[barcode_column]])
  add_df <- ext[matched, cols_to_add, drop = FALSE]
  rownames(add_df) <- rownames(meta)

  if (sum(is.na(matched)) > 0) {
    warning("Not all cells from the Seurat object were found in external metadata.")
  }
  if (length(setdiff(ext[[barcode_column]], meta$Barcode)) > 0) {
    warning("Some cells in external metadata were not found in the Seurat object.")
  }

  object <- AddMetaData(object, metadata = add_df)

  plots <- list()
  if (!is.null(col_to_viz) && length(col_to_viz) == 1 && col_to_viz %in% colnames(object[[]])) {
    if ("umap" %in% Reductions(object)) {
      plots$umap <- DimPlot(object, group.by = col_to_viz, label = FALSE, reduction = "umap") +
        ggtitle(paste(col_to_viz, "(UMAP)"))
    }
    if ("tsne" %in% Reductions(object)) {
      plots$tsne <- DimPlot(object, group.by = col_to_viz, label = FALSE, reduction = "tsne") +
        ggtitle(paste(col_to_viz, "(tSNE)"))
    }
  }

  return(list(object = object, plots = plots))
}
