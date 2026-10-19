
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
  input$object@meta.data[["SCT_snn_res.3.0"]] <-
    seq_len(nrow(input$object@meta.data))
  input$old.columns.to.save <-
    c(
      input$old.columns.to.save,
      "SCT_snn_res.3.0"
    )
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
    output$object@meta.data[["old_SCT_snn_res.3.0"]],
    preserved.metadata[["old_SCT_snn_res.3.0"]]
  )
  }
)


test_that("Recluster requires exact metadata column names.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  possible.metadata.columns <- unique(c(
    intersect(
      c("seurat_clusters", "seurat.clusters"),
      colnames(input$object@meta.data)
    ),
    grep("^SCT_snn_res", colnames(input$object@meta.data), value = TRUE)
  ))
  input$old.columns.to.save <- "SCT_snn_res_2_4"
  expected.message <- paste0(
    "Old column(s) to save not found in metadata: 'SCT_snn_res_2_4'.\n",
    "Possible clustering columns to save before replacement:\n",
    paste0("  - ", possible.metadata.columns, collapse = "\n"),
    "\nIf you want to retain old cluster columns, select their exact names ",
    "with old.columns.to.save and use prepend.txt to give them a unique prefix."
  )

  expect_error(
    do.call(reclusterSeuratObject, input),
    expected.message,
    fixed = TRUE
  )
})

test_that("Recluster reports rejected columns without repeating candidates.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  metadata.names <- colnames(input$object@meta.data)
  input$old.columns.to.save <- c(
    "SCT_snn_res.2.4",
    "missing.column_1",
    "missing.column.2",
    "orig.ident"
  )

  expect_message(
    expect_warning(
      output <- do.call(reclusterSeuratObject, input),
      "Old column(s) to save not found in metadata: 'missing.column_1', 'missing.column.2'.",
      fixed = TRUE
    ),
    paste0(
      "Old column(s) to save are not clustering columns replaced during ",
      "reclustering: 'orig.ident'.\n",
      "     These metadata columns are retained automatically."
    ),
    fixed = TRUE
  )
  expect_identical(colnames(input$object@meta.data), metadata.names)
  expect_true("old_SCT_snn_res.2.4" %in% colnames(output$object@meta.data))
})

test_that("Recluster distinguishes invalid old-column selections.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  metadata.names <- colnames(input$object@meta.data)
  possible.metadata.columns <- unique(c(
    intersect(
      c("seurat_clusters", "seurat.clusters"),
      metadata.names
    ),
    grep("^SCT_snn_res", metadata.names, value = TRUE)
  ))
  input$old.columns.to.save <- c("orig.ident", "missing.column")

  expect_error(
    do.call(reclusterSeuratObject, input),
    paste0(
      "Old column(s) to save not found in metadata: 'missing.column'.\n",
      "Old column(s) to save are not clustering columns replaced during ",
      "reclustering: 'orig.ident'.\n",
      "     These metadata columns are retained automatically.\n",
      "Possible clustering columns to save before replacement:\n",
      paste0("  - ", possible.metadata.columns, collapse = "\n"),
      "\nIf you want to retain old cluster columns, select their exact names ",
      "with old.columns.to.save and use prepend.txt to give them a unique prefix."
    ),
    fixed = TRUE
  )
})


test_that("Recluster supports empty metadata preservation selections.", {

  for (old.columns.to.save in list(NULL, character(0))) {
    input <- getParamsReclusterSeuratObject("Chariou")
    input$object@meta.data[["seurat.clusters"]] <- "legacy"
    input$object@meta.data[["seurat_clusters"]] <- "canonical"
    old1.cluster.values <- rep(c("0", "1"), length.out = nrow(input$object@meta.data))
    input$object@meta.data[["SCT_snn_res.0.2"]] <- old1.cluster.values
    input$object@meta.data[["old1_SCT_snn_res.0.2"]] <- old1.cluster.values
    possible.metadata.columns <- unique(c(
      intersect(
        c("seurat_clusters", "seurat.clusters"),
        colnames(input$object@meta.data)
      ),
      grep("^SCT_snn_res", colnames(input$object@meta.data), value = TRUE)
    ))
    input$number.of.pcs <- if (is.null(old.columns.to.save)) 0 else 10
    input$cluster.resolution.low.range <- 0.4
    input$cluster.resolution.high.range <- 0.4
    input$cluster.resolution.range.bins <- 0.1

    input$old.columns.to.save <- NULL
    output <- NULL
    messages <- character()
    withCallingHandlers(
      output <- do.call(
        reclusterSeuratObject,
        c(input, list(old.columns.to.save = old.columns.to.save))
      ),
      message = function(condition) {
        messages <<- c(messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      }
    )
    expect_true("No old columns were saved.\n" %in% messages)
    expect_false(any(grepl("Possible clustering columns", messages)))

    expect_false("seurat.clusters" %in% colnames(output$object@meta.data))
    expect_false(any(
      setdiff(possible.metadata.columns, "seurat_clusters") %in%
        colnames(output$object@meta.data)
    ))
    expect_false("SCT_snn_res.0.2" %in% colnames(output$object@meta.data))
    expect_true("old1_SCT_snn_res.0.2" %in% colnames(output$object@meta.data))
    expect_equal(
      output$object@meta.data[["old1_SCT_snn_res.0.2"]],
      old1.cluster.values
    )
    expect_true("SCT_snn_res.0.4" %in% colnames(output$object@meta.data))
    expect_true("seurat_clusters" %in% colnames(output$object@meta.data))
    expect_true("sample_name" %in% colnames(output$object@meta.data))
    expect_false(any(grepl("^old_SCT_snn_res", colnames(output$object@meta.data))))
  }
})


test_that("Recluster preserves prior reclustering columns when saving current columns.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  current.cluster.values <- rep(
    c("0", "1"),
    length.out = nrow(input$object@meta.data)
  )
  old1.cluster.values <- rep(
    c("2", "3"),
    length.out = nrow(input$object@meta.data)
  )
  input$object@meta.data[["SCT_snn_res.0.2"]] <- current.cluster.values
  input$object@meta.data[["old1_SCT_snn_res.0.2"]] <- old1.cluster.values
  input$old.columns.to.save <- "SCT_snn_res.0.2"
  input$prepend.txt <- "old"
  input$number.of.pcs <- 10
  input$cluster.resolution.low.range <- 0.4
  input$cluster.resolution.high.range <- 0.4
  input$cluster.resolution.range.bins <- 0.1

  output <- do.call(reclusterSeuratObject, input)

  expect_false("SCT_snn_res.0.2" %in% colnames(output$object@meta.data))
  expect_equal(
    output$object@meta.data[["old_SCT_snn_res.0.2"]],
    current.cluster.values
  )
  expect_equal(
    output$object@meta.data[["old1_SCT_snn_res.0.2"]],
    old1.cluster.values
  )
})


test_that("Recluster rejects existing renamed clustering column names.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  input$object@meta.data[["SCT_snn_res.0.2"]] <- "current"
  input$object@meta.data[["existing_SCT_snn_res.0.2"]] <- "prior"
  input$old.columns.to.save <- "SCT_snn_res.0.2"
  input$prepend.txt <- "existing"

  expect_error(
    do.call(reclusterSeuratObject, input),
    paste0(
      "metadata column names already exist: 'existing_SCT_snn_res.0.2'. ",
      "Set prepend.txt to an unused prefix so the renamed ",
      "column name is unique."
    ),
    fixed = TRUE
  )
})


test_that("Recluster supports custom prefixes, PCA plots, and resolutions.", {

  input <- getParamsReclusterSeuratObject("Chariou")
  input$object <- FindNeighbors(
    input$object,
    reduction = "pca",
    dims = 1:10,
    graph.name = c("harmony_nn", "harmony_snn"),
    verbose = FALSE
  )
  input$object@meta.data[["harmony_snn_res.2.4"]] <-
    input$object@meta.data[["SCT_snn_res.2.4"]]
  input$object@meta.data[["SCT_snn_res_2_4"]] <-
    input$object@meta.data[["SCT_snn_res.2.4"]]
  input$object@meta.data[["SCT-snn-res:2.4"]] <-
    input$object@meta.data[["SCT_snn_res.2.4"]]
  input$old.columns.to.save <- c(
    "SCT_snn_res.2.4",
    "harmony_snn_res.2.4",
    "SCT_snn_res_2_4",
    "SCT-snn-res:2.4"
  )
  input$prepend.txt <- "prior.batch-1.2"
  input$number.of.pcs <- 10
  input$reduction.type <- "pca"
  input$cluster.resolution.low.range <- 0.4
  input$cluster.resolution.high.range <- 0.4
  input$cluster.resolution.range.bins <- 0.1
  expected.sct.values <- input$object@meta.data[["SCT_snn_res.2.4"]]
  expected.harmony.values <- input$object@meta.data[["harmony_snn_res.2.4"]]
  expected.underscore.values <- input$object@meta.data[["SCT_snn_res_2_4"]]
  expected.mixed.separator.values <- input$object@meta.data[["SCT-snn-res:2.4"]]

  messages <- character()
  output <- withCallingHandlers(
    do.call(reclusterSeuratObject, input),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )
  cluster.columns <- grep("^SCT_snn_res\\.", colnames(output$object@meta.data), value = TRUE)

  expect_true(any(grepl(
    "SCT_snn_res.2.4 -> prior.batch-1.2_SCT_snn_res.2.4",
    messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "harmony_snn_res.2.4 -> prior.batch-1.2_harmony_snn_res.2.4",
    messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "SCT_snn_res_2_4 -> prior.batch-1.2_SCT_snn_res_2_4",
    messages,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "SCT-snn-res:2.4 -> prior.batch-1.2_SCT-snn-res:2.4",
    messages,
    fixed = TRUE
  )))
  removal.messages <- messages[grepl(
    "Removing original clustering columns",
    messages,
    fixed = TRUE
  )]
  expect_false(any(grepl(
    "  - SCT_snn_res.2.4",
    removal.messages,
    fixed = TRUE
  )))
  expect_equal(
    output$object@meta.data[["prior.batch-1.2_SCT_snn_res.2.4"]],
    expected.sct.values
  )
  expect_equal(
    output$object@meta.data[["prior.batch-1.2_harmony_snn_res.2.4"]],
    expected.harmony.values
  )
  expect_equal(
    output$object@meta.data[["prior.batch-1.2_SCT_snn_res_2_4"]],
    expected.underscore.values
  )
  expect_equal(
    output$object@meta.data[["prior.batch-1.2_SCT-snn-res:2.4"]],
    expected.mixed.separator.values
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
