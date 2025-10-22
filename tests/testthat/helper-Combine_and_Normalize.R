getParamCN <- function(data){
  
  if(data == "TEC"){
    object <- readRDS(test_path(
      paste0("fixtures/",data), paste0(data,'_Filtered_SO_downsample.rds')))
    object=object
    
    npcs = 30
    vars.to.regress = NULL
    clust.res.low=0.2
    clust.res.high = 1.2
    only.var.genes = FALSE 
    
    
  } else if (data == "Chariou") {
    object <- readRDS(test_path(
      paste0("fixtures/",data), paste0(data,'_Filtered_SO_downsample.rds')))
    object=object
  
    npcs = 21
    vars.to.regress = c("percent.mt")
    clust.res.low=2.4
    clust.res.high = 2.8
    only.var.genes = FALSE 
    
    
  } else if (data == "NSCLC_Single") {
    object <- readRDS(test_path(
      paste0("fixtures/",data), paste0('NSCLCsingle','_Filtered_SO_downsample.rds')))
    
    npcs = 30
    vars.to.regress = NULL
    clust.res.low=0.2
    clust.res.high = 1.2
    only.var.genes = FALSE 
    
    
  } else if (data == "NSCLC_Multi") {
    object <- readRDS(test_path(
      paste0("fixtures/",data), paste0('NSCLC_multi','_Filtered_SO_downsample.rds')))
    object=object
    
    npcs = 30
    vars.to.regress = NULL
    clust.res.low=0.2
    clust.res.high = 1.2
    only.var.genes = FALSE 
    
    
  } else if (data == "BRCA") {
    # object <- readRDS(test_path(
    #   paste0("fixtures/",data), paste0('NSCLCmulti','_Filtered_SO_downsample.rds')))
    # object=object
    # 
    # npcs = 30
    # vars.to.regress = NULL
    # clust.res.low=0.2
    # clust.res.high = 1.2
    # only.var.genes = TRUE 
    
  }
  
  return(list("object" = object,
              "npcs" = npcs,
              "vars.to.regress" = vars.to.regress,
              "clust.res.low"=clust.res.low,
              "clust.res.high" = clust.res.high,
              "only.var.genes" = only.var.genes ))  
}



.drawFig <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".png")
  ggsave(path, x, width = 10, height = 10)
  print(path)
}
.saveSO <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".rds")
  saveRDS(x, file = path)
  print(path)
}

