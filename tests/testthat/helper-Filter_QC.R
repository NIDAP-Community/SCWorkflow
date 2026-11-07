getParamFQ <- function(data){
  
  if(data == "TEC"){
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,'_ProcessRaw_SO_downsample.rds')))
    mad.mitoch.limits=c(NA,3)
    mitoch.limits = c(NA,10)

    
    
  } else if (data == "Chariou") {
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,'_ProcessRaw_SO_downsample.rds')))
    mad.mitoch.limits=c(NA,3)
    mitoch.limits = c(NA,25)

    
    
  } else if (data == "NSCLC_Single") {
    filerds=test_path(
      "fixtures/NSCLC_Single", "NSCLC_Single_ProcessRaw_SO_downsample.rds")
    print(filerds)
    object=readRDS(filerds)

    mad.mitoch.limits=c(NA,3)
    mitoch.limits = c(NA,25)

    
    
  } else if (data == "NSCLC_Multi") {
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,'_ProcessRaw_SO_downsample.rds')))

    mad.mitoch.limits=c(NA,3)
    mitoch.limits = c(NA,25)
    
    
  } else if (data == "BRCA") {
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,'_ProcessRaw_SO_downsample.rds')))
    
    mad.mitoch.limits=c(NA,3)
    mitoch.limits = c(NA,25)
    
    
  }
  
  return(list("object" = object, 
              "mad.mitoch.limits" = mad.mitoch.limits,
              "mitoch.limits" = mitoch.limits,
              "do.doublets.filter" = TRUE))  
}
