getParam3D <- function(data) {
  if (data == "TEC") {
    object <- selectCRObject("TEC")
    color.variable = "orig.ident"
    label.variable = "SCT_snn_res.0.2"
    
  } else if (data == "Chariou") {
    object <- selectCRObject("Chariou")
    color.variable = "orig.ident"
    label.variable = "SCT_snn_res.2.4"
    
  } else if (data == "NSCLC_Single") {
    object <- selectCRObject("NSCLC_Single")
    color.variable = "SCT_snn_res.0.2"
    label.variable = "orig.ident"
    
  } else if (data == "nsclc-multi") {
    object <- selectCRObject("nsclc-multi")
    color.variable = "orig.ident"
    label.variable = "SCT_snn_res.0.2"
    
  } else if (data == "BRCA") {
    object <- selectCRObject("BRCA")
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