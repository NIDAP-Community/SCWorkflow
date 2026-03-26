
getParamsNameClus <- function(data){

  if(data == "TEC"){
    object <- selectSRObject("TEC")
    cluster.numbers <- c(0,1,2,3,4,5,6,7,8,9)
    cluster.names <- c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10")
    cluster.column <- "SCT_snn_res.0.2"
    labels.column <- "mouseRNAseq_main"
  } else if(data == "Chariou"){
    object <- selectSRObject("Chariou")
    cluster.numbers <- c(0,1,2,3,4,5,6,7,8,9,10,
                         11,12,13,14,15,16,17,18,19,20,
                         21,22,23,24,25,26,27,28,29,30,
                         31,32,33,34,35,36)
    cluster.names <- c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10",
                       "L11","L12","L13","L14","L15","L16","L17","L18","L19",
                       "L20","L21","L22","L23","L24","L25","L26","L27","L28",
                       "L29","L30","L31","L32","L33","L34","L35","L36","L37")
    cluster.column <- "SCT_snn_res.2.4"
    labels.column <- "immgen_main"
  } else if(data == "pbmc-single"){
      object <- selectSRObject("pbmc-single")
      cluster.numbers <- c(0,1,2,3,4,5,6,7,8,9,10,
                           11,12,13)
      cluster.names <- c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10",
                         "L11","L12","L13","L14")
      cluster.column <- "SCT_snn_res.0.2"
      labels.column <- "BP_encode_main"
  } else if(data == "nsclc-multi"){
      object <- selectSRObject("nsclc-multi")
      cluster.numbers <- c(0,1,2,3,4,5,6,7,8,9,10,
                           11,12,13,14,15)
      cluster.names <- c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10",
                         "L11","L12","L13","L14","L15","L16")
      cluster.column <- "SCT_snn_res.0.2"
      labels.column <- "BP_encode"
  }
  
  cluster.identities.table <- data.frame(
    OriginalClusterIDs = cluster.numbers,
    NewClusterNames = cluster.names,
    stringsAsFactors = FALSE
  )

  return(list("object" = object,
              "cluster.column" = cluster.column,
              "labels.column" = labels.column,
              "cluster.identities.table" = cluster.identities.table,
              "cluster.numbers" = "OriginalClusterIDs",
              "cluster.names" = "NewClusterNames"
              ))
  
}

.drawpng <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".png")
  png(path,
      width=width,
      height=height,
      units = "in",
      res = 400
  )
  on.exit(dev.off())
  print(x)
  path
}

