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

.drawplot <- function(x) {
  path <- tempfile(fileext = ".png")
  plotly::save_image(x, path)
  path
}