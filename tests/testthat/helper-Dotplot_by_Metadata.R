getParamDP <- function(data){

  if(data == "TEC"){
      object <- selectCRObject("TEC")
      metadata <- "SCT_snn_res.0.2"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c(0,1,2,3,4,5,6,7,8,9)
  } else if (data == "Chariou") {
      object <- selectCRObject("Chariou")
      metadata <- "orig.ident"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("PBS","CD8dep","NHSIL12","ENT","Combo")
  } else if (data == "NSCLC_Single") {
      object <- selectSRObject("NSCLC_Single")
      metadata <- "BP_encode_main"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("B-cells","CD4+ T-cells","CD8+ T-cells","DC",
                 "HSC","Monocytes","NK cells")  
  } else if (data == "nsclc-multi") {
      object <- selectSRObject("nsclc-multi")
      metadata <- "BP_encode_main"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("B-cells","CD4+ T-cells","CD8+ T-cells","DC",
                 "HSC","Monocytes","NK cells")  
  } else if (data == "BRCA") {
      object <- selectCRObject("BRCA")
      metadata <- "SCT_snn_res.0.2"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("0","1","2","3","4","5","6")  
  }
  
  return(list("object" = object, 
              "metadata"=metadata, 
              "markers" = markers, 
              "cells" = cells)
         )  
}

.drawDotPng <- function(x, width = 10, height = 10){
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