# Helper functions for Compare_Cell_Populations tests

# Load real Seurat objects from fixtures
getParamCCP <- function(data) {
  if (data == "TEC") {
    object <- selectCRObject("TEC")
    annotation.column <- "seurat_clusters"
    group.column <- "Status"
    sample.column <- "orig.ident"
    counts.type <- "Frequency"
    group.order <- NULL
    
  } else if (data == "Chariou") {
    object <- selectCRObject("Chariou")
    annotation.column <- "seurat_clusters"
    group.column <- "Status"
    sample.column <- "orig.ident"
    counts.type <- "Frequency"
    group.order <- NULL
    
  } else if (data == "PBMC") {
    object <- selectSRObject("pbmc-single")
    annotation.column <- "HPCA_main"
    group.column <- "Phase"
    sample.column <- "orig.ident"
    counts.type <- "Frequency"
    group.order <- NULL
    
  } else if (data == "NSCLC") {
    object <- selectCRObject("nsclc-multi")
    annotation.column <- "seurat_clusters"
    group.column <- "Phase"
    sample.column <- "orig.ident"
    counts.type <- "Frequency"
    group.order <- NULL
    
  } else if (data == "BRCA") {
    object <- selectCRObject("BRCA")
    annotation.column <- "seurat_clusters"
    group.column <- "Phase"
    sample.column <- "orig.ident"
    counts.type <- "Frequency"
    group.order <- NULL
  } else {
    stop(
      "Unknown dataset '", data, "'. ",
      "Supported datasets are: 'TEC', 'Chariou', 'PBMC', 'NSCLC', 'BRCA'."
    )
  }
  
  return(
    list(
      "object" = object,
      "metadata.table" = object@meta.data,
      "annotation.column" = annotation.column,
      "group.column" = group.column,
      "sample.column" = sample.column,
      "counts.type" = counts.type,
      "group.order" = group.order
    )
  )
}

# Helper function to save ggplot objects for snapshot testing
.drawCCPFig <- function(x, width = 10, height = 10) {
  path <- tempfile(fileext = ".png")
  ggplot2::ggsave(path, x, width = width, height = height)
  path
}

# Helper function to save data tables for snapshot testing
.saveCCPTable <- function(x) {
  path <- tempfile(fileext = ".rds")
  saveRDS(x, file = path)
  path
}
