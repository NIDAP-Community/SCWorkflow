getCbgAutoParam <- function(data) {
  if (data == "TEC") {
    object = selectCRObject("TEC")
    samples.subset = unique(object$orig.ident)
    samples.to.display = unique(object$orig.ident)
    marker.table = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
    cells.of.interest = colnames(marker.table)[1:3]
    
  } else if (data == "Chariou") {
    object = selectCRObject("Chariou")
    samples.subset = unique(object$orig.ident)
    samples.to.display = unique(object$orig.ident)
    marker.table = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
    cells.of.interest = colnames(marker.table)[4:6]
    
  } else if (data == "pbmc.single") {
    object = selectCRObject("pbmc-single")
    samples.subset = unique(object$orig.ident)
    samples.to.display = unique(object$orig.ident)
    marker.table = apply(read.csv(test_path("fixtures", 
                                            "Marker_Table_demo.csv")), 
                         2, toupper)
    set.seed(93)
    marker.table = data.frame(
      rand_type1 = sample(rownames(object), 5, replace = FALSE),
      rand_type2 = sample(rownames(object), 5, replace = FALSE),
      rand_type3 = sample(rownames(object), 5, replace = FALSE)
    )
    cells.of.interest = colnames(marker.table)
    
  } else if (data == "nsclc.multi") {
    object = selectCRObject("nsclc-multi")
    samples.subset = unique(object$orig.ident)
    samples.to.display = unique(object$orig.ident)
    set.seed(94)
    marker.table = data.frame(
      rand_type1 = sample(rownames(object), 5, replace = FALSE),
      rand_type2 = sample(rownames(object), 5, replace = FALSE),
      rand_type3 = sample(rownames(object), 5, replace = FALSE)
    )
    cells.of.interest = colnames(marker.table)
    
  } else if (data == "BRCA") {
    object = selectCRObject("BRCA")
    samples.subset = unique(object$orig.ident)
    samples.to.display = unique(object$orig.ident)
    set.seed(95)
    marker.table = data.frame(
      rand_type1 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE),
      rand_type2 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE),
      rand_type3 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE)
    )
    cells.of.interest = colnames(marker.table)
    
  }
  
  return(
    list(
      "object" = object,
      "samples.subset" = samples.subset,
      "samples.to.display" = samples.to.display,
      "marker.table" = marker.table,
      "cells.of.interest" = cells.of.interest
    )
  )
}

.drawCbG <- function(x, width = 10, height = 10, component = c("overall", "celltype", "manual_entry"), index = 1){
  component <- match.arg(component)
  target <- x
  if (is.list(x) && all(c("overall", "celltype", "manual_entry") %in% names(x))) {
    target <- x[[component]]
    if (is.list(target)) {
      if (length(target) < index) {
        stop("Requested index exceeds available plots in component")
      }
      target <- target[[index]]
    }
    if (is.null(target)) {
      # Fallback to first available overall plot if chosen component is NULL
      if (length(x$overall) >= 1) {
        target <- x$overall[[1]]
      } else {
        stop("No plot available to save in provided object")
      }
    }
  }
  path <- tempfile(fileext = ".png")
  ggsave(path, target, width = width, height = height)
  print(path)
}
