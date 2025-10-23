selectCRObject <- function(dataset) {
  
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    input.object <- readRDS(test_path("fixtures/TEC", 
        "TEC_Combine_and_Renormalize_SO_downsample.rds"))
    
  } else if (dataset == "Chariou"){
    
    print("selected Chariou dataset")
    input.object <- readRDS(test_path("fixtures/Chariou", 
        "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
    
  } else if (dataset == "pbmc-single"){
    
    print("selected PBMC Single dataset")
    input.object <- readRDS(test_path("fixtures/PBMC_Single",
    "PBMCsingle_Cell_Types_SingleR_SO_downsample.rds"))
    
  } else if (dataset == "nsclc-multi"){
    
    print("selected NSCLC Multi dataset")
    input.object <- readRDS(test_path("fixtures/NSCLC_Multi",
    "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
  
  } else if (dataset == "BRCA"){
    
    print("selected BRCA dataset")
    input.object <- readRDS(test_path("fixtures/BRCA",
        "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
  }
  return(input.object) 
  
}


selectSRObject <- function(dataset){
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    input.object <- readRDS(test_path("fixtures/TEC", 
          "TEC_CellTypesSingleR_SO_downsample.rds"))
    
  } else if (dataset == "Chariou"){
    
    print("selected Chariou dataset")
    input.object <- readRDS(test_path("fixtures/Chariou", 
          "Chariou_Cell_Types_SingleR_SO_downsample.rds"))
    
  } else if (dataset == "pbmc-single"){
    
    print("selected PBMCsingle dataset")
    input.object <- readRDS(test_path("fixtures/PBMC_Single", 
          "PBMCsingle_Cell_Types_SingleR_SO_downsample.rds"))
  
  } else if (dataset == "nsclc-multi"){
    
    print("selected NSCLCmulti dataset")
    input.object <- readRDS(test_path("fixtures/NSCLC_Multi", 
       "NSCLCmulti_Cell_Types_SingleR_SO_downsample.rds"))
  
  } else if (dataset == "BRCA"){
    
    print("selected BRCA dataset")
    input.object <- readRDS(test_path("fixtures/BRCA", 
                                      "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
    
  }
}
  
