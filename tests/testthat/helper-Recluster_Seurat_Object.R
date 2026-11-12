
getParamsReclusterSeuratObject <- function(data) {
  
  ## If TEC (Mouse) dataset selected:
  if (data == "TEC") {
    object <- selectSRObject("TEC")
    old.columns.to.save = c(
      "SCT_snn_res.0.2",
      "SCT_snn_res.0.4",
      "SCT_snn_res.0.6",
      "SCT_snn_res.0.8",
      "SCT_snn_res.1",
      "SCT_snn_res.1.2"
    )

  ## Else if Charoiu (Mouse) dataset selected:
  } else if (data == "Chariou") {
    object <- selectSRObject("Chariou")
    old.columns.to.save = c(
      "SCT_snn_res.2.4",
      "SCT_snn_res.2.6",
      "SCT_snn_res.2.8"
    )
    
  ## Else if NSCLC_Single dataset selected:
  } else if (data == "NSCLC_Single") {
    object <- selectSRObject("NSCLC_Single")
    old.columns.to.save = c(
      "SCT_snn_res.0.2",
      "SCT_snn_res.0.4",
      "SCT_snn_res.0.6",
      "SCT_snn_res.0.8",
      "SCT_snn_res.1",
      "SCT_snn_res.1.2"
    )
  
  ## Else if NSCLCMulti dataset selected:
  } else if (data == "nsclc-multi") {
    object <- selectSRObject("nsclc-multi")
    old.columns.to.save = c(
      "SCT_snn_res.0.2",
      "SCT_snn_res.0.4",
      "SCT_snn_res.0.6",
      "SCT_snn_res.0.8",
      "SCT_snn_res.1",
      "SCT_snn_res.1.2"
    )
  
  ## Else if BRCA dataset selected:
  } else if (data == "BRCA") {
    object <- selectCRObject("BRCA")
    old.columns.to.save = c(
      "SCT_snn_res.0.2",
      "SCT_snn_res.0.4",
      "SCT_snn_res.0.6",
      "SCT_snn_res.0.8",
      "SCT_snn_res.1",
      "SCT_snn_res.1.2"
    )
  }

return(list("object" = object,
            "old.columns.to.save" = old.columns.to.save))
}
