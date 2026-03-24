
test_that("Test Recluster with TEC (Mouse) dataset.", {
  
  ## Get inputs and default parameters.
  input <- getParamsReclusterSeuratObject("TEC")
  
  ## Do function call.
  output <- do.call(reclusterSeuratObject, input)
  
  ## Compare new plot to old plot.
  ggsave(
    "output/TEC_reclusterSO.png",
    output$plots,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_reclusterSO.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)


test_that("Test Recluster with TEC (Mouse) dataset & UMAP.", {
  
  ## Get inputs and default parameters.
  input <- getParamsReclusterSeuratObject("TEC")
  
  ## Set non-default parameters for this test.
  input$reduction.type = "umap"
  
  ## Do function call.
  output <- do.call(reclusterSeuratObject, input)
  
  ## Compare new plot to old plot.
  ggsave(
    "output/TEC_reclusterSO_umap.png",
    output$plots,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_reclusterSO_umap.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)


test_that("Test Recluster with Chariou (Mouse) dataset.", {

  ## Get inputs and default parameters.
  input <- getParamsReclusterSeuratObject("Chariou")
  
  ## Do function call.
  output <- do.call(reclusterSeuratObject, input)
  
  ## Compare new plot to old plot.
  ggsave(
    "output/Chariou_reclusterSO.png",
    output$plots,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "Chariou_reclusterSO.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)


test_that("Test Recluster with PBMCSingle (Human) dataset.", {
  
  ## Get inputs and default parameters.
  input <- getParamsReclusterSeuratObject("pbmc-single")
  
  ## Do function call.
  output <- do.call(reclusterSeuratObject, input)
  
  ## Compare new plot to old plot.
  ggsave(
    "output/PBMCsingle_reclusterSO.png",
    output$plots,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "PBMCsingle_reclusterSO.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)


test_that("Test Recluster with NSCLCMulti (Human) dataset.", {
  
  ## Get inputs and default parameters.
  input <- getParamsReclusterSeuratObject("nsclc-multi")
  
  ## Do function call.
  output <- do.call(reclusterSeuratObject, input)
  
  ## Compare new plot to old plot.
  ggsave(
    "output/NSCLCMulti_reclusterSO.png",
    output$plots,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "NSCLCMulti_reclusterSO.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)


test_that("Test Recluster with BRCA (Human) dataset.", {
  
  ## Get inputs and default parameters.
  input <- getParamsReclusterSeuratObject("BRCA")
  
  ## Do function call.
  output <- do.call(reclusterSeuratObject, input)
  
  ## Compare new plot to old plot.
  ggsave(
    "output/BRCA_reclusterSO.png",
    output$plots,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "BRCA_reclusterSO.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)
