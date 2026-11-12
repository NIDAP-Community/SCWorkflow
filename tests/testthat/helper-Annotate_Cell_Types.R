#This helper script will return parameters for each test dependent on the input:

getParamACT <- function(data){

  if(data == "TEC"){
      object <- selectCRObject("TEC")
      species <- "Mouse"
      use.clusters = "seurat_clusters"
  } else if (data == "Chariou") {
      object <- selectCRObject("Chariou")
      species <- "Mouse"
      use.clusters = "seurat_clusters"
  } else if (data == "NSCLC_Single") {
      object <- selectCRObject("NSCLC_Single")
      species <- "Human"
      use.clusters <- NULL
  } else if (data == "nsclc-multi") {
      object <- selectCRObject("nsclc-multi")
      species <- "Human"
      use.clusters = "seurat_clusters"
  } else if (data == "BRCA") {
      object <- selectCRObject("BRCA")
      species <- "Human"
      use.clusters = "seurat_clusters"
  }
 
  return(list("object" = object, "species" = species, "use.clusters" = use.clusters))  
}
