
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
  expect_snapshot_file("output/TEC_reclusterSO.png")
  
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
  expect_snapshot_file("output/TEC_reclusterSO_umap.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)


test_that("Test Recluster with Chariou (Mouse) dataset.", {

  ## Get inputs and default parameters.
  input <- getParamsReclusterSeuratObject("Chariou")
  preserved.metadata <- input$object@meta.data[, input$old.columns.to.save, drop = FALSE]
  preserved.names <- paste0("old_", input$old.columns.to.save)
  colnames(preserved.metadata) <- preserved.names
  
  ## Do function call.
  output <- do.call(reclusterSeuratObject, input)
  
  ## Compare new plot to old plot.
  ggsave(
    "output/Chariou_reclusterSO.png",
    output$plots,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output/Chariou_reclusterSO.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  expect_true(all(preserved.names %in% colnames(output$object@meta.data)))
  expect_equal(output$object@meta.data[, preserved.names, drop = FALSE], preserved.metadata)
  }
)


test_that("Recluster requires exact metadata column names.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  input$old.columns.to.save <- "SCT_snn_res_2_4"

  expect_error(
    do.call(reclusterSeuratObject, input),
    "Could not find requested metadata columns: SCT_snn_res_2_4."
  )
})


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
  expect_snapshot_file("output/PBMCsingle_reclusterSO.png")
  
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
  expect_snapshot_file("output/NSCLCMulti_reclusterSO.png")
  
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
  expect_snapshot_file("output/BRCA_reclusterSO.png")
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)
