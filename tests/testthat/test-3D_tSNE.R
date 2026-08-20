#These lines are necessary to launch Kaleido properly (temporary patch):
skip_on_ci()
reticulate::py_require(c("plotly", "kaleido==0.2.1"))
reticulate::py_run_string(
  "import sys;print(sys.version); sys.path.append('/rstudio-files/R_environments/single-cell-rna-seq-r4'); print(sys.path)"
)

test_that("Produce 3D tsne plot and return tsne coordinates - TEC Data", {
  cr.object <- getParam3D("TEC")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)
  
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)
 
  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plots),
    "TEC_plotly.png"
    )
  }
)

test_that("Run 3DTSNE with error for color selection - TEC Data", {
  cr.object <- getParam3D("TEC")
  cr.object$plot.type <- "tSNE"
  cr.object$color.variable <- "Likely_CellType"
  expect_error(output <- do.call(tSNE3D, cr.object),
               "^The metadata variable selected for color")

})

test_that("Run 3DTSNE with error for color selection - TEC Data", {
  cr.object <- getParam3D("TEC")
  cr.object$plot.type <- "tSNE"
  cr.object$label.variable <- "Likely_CellType"
  expect_error(output <- do.call(tSNE3D, cr.object),
               "^The metadata variable selected for labeling")

})

test_that("Produce 3D tsne plot and return tsne coordinates - Chariou Data", {
  cr.object <- getParam3D("Chariou")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)

  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)
  expect_true(all(c("tSNE_1", "tSNE_2", "tSNE_3") %in% colnames(output$data)))

  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plots),
    "Chariou_plotly.png"
  )
}
)

test_that("Produce 3D UMAP plot and return UMAP coordinates - Chariou Data", {
  cr.object <- getParam3D("Chariou")
  cr.object$plot.type <- "UMAP"
  output <- do.call(tSNE3D, cr.object)

  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)
  expect_true(all(c("UMAP_1", "UMAP_2", "UMAP_3") %in% colnames(output$data)))
})

test_that("Run 3D plot with error for invalid plot type", {
  cr.object <- getParam3D("Chariou")
  cr.object$plot.type <- "PCA"
  expect_error(do.call(tSNE3D, cr.object), "plot.type must be one of")
})

test_that("Produce 3D tsne plot and return tsne - PBMC-single Data", {
  cr.object <- getParam3D("pbmc-single")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plots),
    "PBMC_single_plotly.png"
  )
})

test_that("Produce 3D tsne plot and return tsne - NSCLC-multi Data", {
  cr.object <- getParam3D("nsclc-multi")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)

  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plots),
    "NSCLC_multi_plotly.png"
  )
})

test_that("Produce 3D tsne plot and return tsne - BRCA Data", {
  cr.object <- getParam3D("BRCA")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)

  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plots),
    "BRCA_plotly.png"
  )
})
