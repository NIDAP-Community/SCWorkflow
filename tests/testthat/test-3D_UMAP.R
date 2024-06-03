#These lines are necessary to launch Kaleido properly (temporary patch):
skip_on_ci()
reticulate::py_run_string(
  "import sys;print(sys.version); sys.path.append('/rstudio-files/R_environments/single-cell-rna-seq-r4'); print(sys.path)"
)

test_that("Produce 3D UMAP plot and return tsne coordinates - TEC Data", {
  cr.object <- getParam3D("TEC")
  output <- do.call(UMAP3D, cr.object)
  
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
 
  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plot),
    "TEC_plotly.png"
    )
  }
)

test_that("Run 3DUMAP with error for color selection - TEC Data", {
  cr.object <- getParam3D("TEC")
  cr.object$color.variable <- "Likely_CellType"
  expect_error(output <- do.call(UMAP3D, cr.object),
               "^The metadata variable selected for color")

})

test_that("Run 3DUMAP with error for color selection - TEC Data", {
  cr.object <- getParam3D("TEC")
  cr.object$label.variable <- "Likely_CellType"
  expect_error(output <- do.call(UMAP3D, cr.object),
               "^The metadata variable selected for labeling")

})

test_that("Produce 3D UMAP plot and return tsne coordinates - Chariou Data", {
  cr.object <- getParam3D("Chariou")
  output <- do.call(UMAP3D, cr.object)

  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plot),
    "Chariou_plotly.png"
  )
}
)

test_that("Produce 3D UMAP plot and return tsne - PBMC-single Data", {
  cr.object <- getParam3D("pbmc-single")
  output <- do.call(UMAP3D, cr.object)
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plot),
    "PBMC_single_plotly.png"
  )
})

test_that("Produce 3D UMAP plot and return tsne - NSCLC-multi Data", {
  cr.object <- getParam3D("nsclc-multi")
  output <- do.call(UMAP3D, cr.object)

  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plot),
    "NSCLC_multi_plotly.png"
  )
})

test_that("Produce 3D UMAP plot and return tsne - BRCA Data", {
  cr.object <- getParam3D("BRCA")
  output <- do.call(UMAP3D, cr.object)

  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawplot(output$plot),
    "BRCA_plotly.png"
  )
})
