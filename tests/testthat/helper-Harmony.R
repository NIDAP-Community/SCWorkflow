getHarmonyParam <- function(data){

  if(data == "TEC"){
    
      object = selectCRObject("TEC")
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.var = c("orig.ident")
      nvar = 100

  } else if (data == "Chariou") {
    
      object = selectCRObject("Chariou")
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.var = c("orig.ident")
      nvar = 100
      
  } else if (data == "NSCLC_Single") {
    
      object = selectCRObject("NSCLC_Single")
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.var = c("Phase")
      nvar = 100
    
  } else if (data == "nsclc_multi") {

      object = selectCRObject("nsclc-multi")
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.var = c("orig.ident")
      nvar = 100

  } else if (data == "BRCA") {

      object = selectCRObject("BRCA")
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.var = c("orig.ident")
      nvar = 100
    
  }
  
  return(list("object" = object, 
              "genes.to.add"= genes.to.add, 
              "group.by.var" = group.by.var,
              "nvar" = nvar))  
}

.drawHarmonyFig <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".png")
  ggsave(path, x, width = 10, height = 10)
  print(path)
}
