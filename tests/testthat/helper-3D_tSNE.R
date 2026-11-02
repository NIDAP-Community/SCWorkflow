getParam3D <- function(data) {
  if (data == "TEC") {
    object <- selectSRObject("TEC")
    color.variable = "orig.ident"
    label.variable = "immgen_main"
    
  } else if (data == "Chariou") {
    object <- selectSRObject("Chariou")
    color.variable = "mouseRNAseq"
    label.variable = "Barcode"
    
  } else if (data == "pbmc-single") {
    object <- selectSRObject("pbmc-single")
    color.variable = "SCT_snn_res.0.2"
    label.variable = "BP_encode_main"
    
  } else if (data == "nsclc-multi") {
    object <- selectSRObject("nsclc-multi")
    color.variable = "orig.ident"
    label.variable = "BP_encode_main"
    
  } else if (data == "BRCA") {
    object <- selectSRObject("BRCA")
    color.variable = "orig.ident"
    label.variable = "SCT_snn_res.0.4"
    
  }
  
  return(
    list(
      "object" = object,
      "color.variable" = color.variable,
      "label.variable" = label.variable
    )
  )
}

expect3DHtmlOutput <- function(plot) {
  testthat::expect_s3_class(plot, "plotly")
  testthat::expect_s3_class(plot, "htmlwidget")
  path <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(plotly::as_widget(plot), path, selfcontained = TRUE)
  testthat::expect_true(file.exists(path))
  testthat::expect_gt(file.info(path)$size, 0)
  invisible(path)
}