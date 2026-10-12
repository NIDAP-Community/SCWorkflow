
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
  if (!nzchar(Sys.getenv("CI"))) {
    expect_snapshot_file("output/TEC_reclusterSO.png")
  }
  
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
  if (!nzchar(Sys.getenv("CI"))) {
    expect_snapshot_file("output/TEC_reclusterSO_umap.png")
  }
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)


test_that("Test Recluster with Chariou (Mouse) dataset.", {

  ## Get inputs and default parameters.
  input <- getParamsReclusterSeuratObject("Chariou")
  input$object@meta.data[["batch.name_1.2"]] <-
    seq_len(nrow(input$object@meta.data))
  input$old.columns.to.save <-
    c(input$old.columns.to.save, "batch.name_1.2")
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
  if (!nzchar(Sys.getenv("CI"))) {
    expect_snapshot_file("output/Chariou_reclusterSO.png")
  }
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  expect_true(all(preserved.names %in% colnames(output$object@meta.data)))
  expect_equal(output$object@meta.data[, preserved.names, drop = FALSE], preserved.metadata)
  expect_equal(
    output$object@meta.data[["old_batch.name_1.2"]],
    preserved.metadata[["old_batch.name_1.2"]]
  )
  expect_true("batch.name_1.2" %in% colnames(output$object@meta.data))
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


test_that("Recluster reports all missing metadata columns without mutation.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  metadata.names <- colnames(input$object@meta.data)
  input$old.columns.to.save <- c(
    "SCT_snn_res.2.4",
    "missing.column_1",
    "missing.column.2"
  )

  expect_error(
    do.call(reclusterSeuratObject, input),
    paste0(
      "Could not find requested metadata columns: ",
      "missing.column_1, missing.column.2."
    )
  )
  expect_identical(colnames(input$object@meta.data), metadata.names)
})


test_that("Recluster supports empty metadata preservation selections.", {

  for (old.columns.to.save in list(NULL, character(0))) {
    input <- getParamsReclusterSeuratObject("Chariou")
    input$number.of.pcs <- if (is.null(old.columns.to.save)) 0 else 10
    input$cluster.resolution.low.range <- 0.4
    input$cluster.resolution.high.range <- 0.4
    input$cluster.resolution.range.bins <- 0.1

    input$old.columns.to.save <- NULL
    output <- do.call(
      reclusterSeuratObject,
      c(input, list(old.columns.to.save = old.columns.to.save))
    )

    expect_true("SCT_snn_res.0.4" %in% colnames(output$object@meta.data))
    expect_false(any(grepl("^old_SCT_snn_res", colnames(output$object@meta.data))))
  }
})


test_that("Recluster supports custom prefixes, PCA plots, and resolutions.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  input$old.columns.to.save <- "SCT_snn_res.2.4"
  input$prepend.txt <- "prior.batch-1.2"
  input$number.of.pcs <- 10
  input$reduction.type <- "pca"
  input$cluster.resolution.low.range <- 0.4
  input$cluster.resolution.high.range <- 0.4
  input$cluster.resolution.range.bins <- 0.1
  expected.values <- input$object@meta.data[[input$old.columns.to.save]]

  output <- do.call(reclusterSeuratObject, input)
  cluster.columns <- grep("^SCT_snn_res\\.", colnames(output$object@meta.data), value = TRUE)

  expect_equal(
    output$object@meta.data[["prior.batch-1.2_SCT_snn_res.2.4"]],
    expected.values
  )
  expect_setequal(cluster.columns, "SCT_snn_res.0.4")
})


test_that("Recluster validates PCA and resolution parameters.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  input$number.of.pcs <- -1
  expect_error(
    do.call(reclusterSeuratObject, input),
    "number.of.pcs must be a non-negative whole number."
  )

  input$number.of.pcs <- 1.5
  expect_error(
    do.call(reclusterSeuratObject, input),
    "number.of.pcs must be a non-negative whole number."
  )

  input$number.of.pcs <- 10
  input$cluster.resolution.low.range <- 0.6
  input$cluster.resolution.high.range <- 0.4
  expect_error(
    do.call(reclusterSeuratObject, input),
    "cluster.resolution.low.range must be less than or equal to "
  )

  input$cluster.resolution.low.range <- 0.4
  input$cluster.resolution.range.bins <- 0
  expect_error(
    do.call(reclusterSeuratObject, input),
    "cluster.resolution.range.bins must be greater than zero."
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
  if (!nzchar(Sys.getenv("CI"))) {
    expect_snapshot_file("output/PBMCsingle_reclusterSO.png")
  }
  
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
  if (!nzchar(Sys.getenv("CI"))) {
    expect_snapshot_file("output/NSCLCMulti_reclusterSO.png")
  }
  
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
  if (!nzchar(Sys.getenv("CI"))) {
    expect_snapshot_file("output/BRCA_reclusterSO.png")
  }
  
  ## Compare results.
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements) 
  }
)
