#This helper script will return parameters for each test dependent on the input:

getParamFSOBM <- function(data) {
  if (data == "TEC") {
    object <- selectCRObject("TEC")
    samples.to.include <- 'c("1_Embryo_13_5","3_Newborn")'
    sample.name <- "orig.ident"
    category.to.filter <- "seurat_clusters"
    values.to.filter <- "3"
  } else if (data == "Chariou") {
    object <- selectCRObject("Chariou")
    samples.to.include = 'c("CD8dep","NHSIL12")'
    sample.name = "orig.ident"
    category.to.filter = "seurat_clusters"
    values.to.filter = "2"
  } else if (data == "NSCLC_Single") {
    object <- selectCRObject("NSCLC_Single")
    samples.to.include = paste0('c("', unique(as.character(object$orig.ident)), '")')
    sample.name = "orig.ident"
    category.to.filter = "SCT_snn_res_0_4"
    values.to.filter = "5"
  } else if (data == "nsclc-multi") {
    object <- selectCRObject("nsclc-multi")
    samples.to.include <- 'c("Donor_1","Donor_4")'
    sample.name <- "orig.ident"
    category.to.filter <- "seurat_clusters"
    values.to.filter <- "16"
  } else if (data == "BRCA") {
    object <- selectCRObject("BRCA")
    samples.to.include <- 'c("CID3586","CID3946", "CID4513" )'
    sample.name <- "orig.ident"
    category.to.filter <- "seurat_clusters"
    values.to.filter <- "1"
  }
  
  
  
  
  return(
    list(
      "object" = object,
      "samples.to.include" = samples.to.include,
      "sample.name" = sample.name,
      "category.to.filter" = category.to.filter,
      "values.to.filter" = values.to.filter
    )
  )
}


