
#' @title Append Metadata to Seurat Object.
#' @description This template appends sample metadata from an input
#' table to a Seurat object, creating new metadata columns labeling all cells 
#' from each sample by the new metadata values.
#' @details This template appends sample metadata from an input
#' table to a Seurat object, creating new metadata columns labeling all cells 
#' from each sample by the new metadata values.
#'
#' @param object The input Seurat Object or list of Seurat Objects to which you 
#' wish to add metadata.
#' @param metadata.to.append A table of sample metadata that you want to append 
#' to the already-existing metadata within the input Seurat Object(s).
#' @param sample.name.column The column of the input metadata.to.append table
#' that contains sample names matching the orig.idents in the input object(s).
#'
#' @importFrom Seurat AddMetaData
#'
#' @export
#' @return Function returns a Seurat Object or Objects with additional metadata
#' columns containing the appended metadata now annotated each cell by sample
#' name (orig.ident).

## Origin: This code adapted from v7 of our Potomac released NIDAP template:
## Append Metadata to Seurat Objects [scRNA-Seq][CCBR] 
## (3abab407-68fd-4f72-b87e-9044e0c4bb18): v7


appendMetadataToSeuratObject <- function(object,
                                         metadata.to.append,
                                         sample.name.column
                                         ) {
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  ## Get the columns of input sample metadata table you want.
  rownames(metadata.to.append) <- metadata.to.append[[sample.name.column]]
  metacols <- colnames(metadata.to.append)
  metacols <- metacols[metacols != sample.name.column]
  
  
  ## Check if sample number is same input SOs and metadata table.

  if(class(object)=='list'){
    
    if(identical(sort(names(object)), 
                 sort(metadata.to.append[[sample.name.column]]))==F){
    stop("ERROR: Your metadata sample names table does not have the same number 
       of samples as you have Seurat objects.")
    }
    
    ## Append column(s) of metadata.to.append to existing metadata in object(s).
      for(i in names(object)){
        pname <- i
        for(x in 1:length(metacols)) {
          object[[i]]=AddMetaData(object[[i]],
                                 metadata=metadata.to.append[pname,metacols[x]],
                                 col.name=metacols[x])
        }
      }  
      
      
  } else if (class(object)=='Seurat') {
  
  if(identical(sort(unique(object@meta.data$orig.ident)), 
               sort(metadata.to.append[[sample.name.column]]))==F){
    stop("ERROR: Your metadata sample names table does not have the same number 
         of samples as you have Seurat objects.")
  }
    
    object@meta.data=merge(object@meta.data,metadata.to.append,
                           by.x='orig.ident',
                           by.y=sample.name.column,
                           all.x=T)

    if(sum(metacols%in%colnames(object@meta.data))==length(metacols)){
      stop(
        cat("ERROR: Not all Metadata was added to Seurat Object. Make sure that 
            metadata column names are not named: \n",
            paste(colnames(object@meta.data),collapse = '\n'))
        )
    }
    
    }
  

## Print logs heading the altered SOs.
  # print("resulting metadata tables:\n")
  # print(lapply(object, function(x) head(x@meta.data)))
  
  ## Return object as single-item list.
  result.list <- list("object" = object)
  return(result.list)
}
