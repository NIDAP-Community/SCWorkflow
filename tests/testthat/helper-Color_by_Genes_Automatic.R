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

.drawCbG <- function(gglist, width = 20, height = 3 * length(gglist)) {
  # Combine the list of plots into a single plot
  combined_plot <- do.call(grid.arrange, c(gglist, ncol = 3))
  
  # Save the combined plot to a temporary file
  path <- tempfile(fileext = ".png")
  ggsave(path, combined_plot, width = width, height = height)
  return(path)
}
