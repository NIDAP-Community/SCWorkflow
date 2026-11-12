test_that("Produce 3D tsne plot and return tsne coordinates - TEC Data", {
  cr.object <- getParam3D("TEC")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)
  
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)
 
  expect3DHtmlOutput(output$plots)
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

  expect3DHtmlOutput(output$plots)
}
)

test_that("Produce 3D UMAP plot and return UMAP coordinates - Chariou Data", {
  cr.object <- getParam3D("Chariou")
  cr.object$plot.type <- "UMAP"
  output <- do.call(tSNE3D, cr.object)

  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)
  expect_true(all(c("UMAP_1", "UMAP_2", "UMAP_3") %in% colnames(output$data)))
  expect3DHtmlOutput(output$plots)
})

test_that("Run 3D plot with error for invalid plot type", {
  cr.object <- getParam3D("Chariou")
  cr.object$plot.type <- "PCA"
  expect_error(do.call(tSNE3D, cr.object), "plot.type must be one of")
})

test_that("Produce 3D tsne plot and return tsne - NSCLC_Single Data", {
  cr.object <- getParam3D("NSCLC_Single")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  expect3DHtmlOutput(output$plots)
})

test_that("Produce 3D tsne plot and return tsne - NSCLC-multi Data", {
  cr.object <- getParam3D("nsclc-multi")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)

  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  expect3DHtmlOutput(output$plots)
})

test_that("Produce 3D tsne plot and return tsne - BRCA Data", {
  cr.object <- getParam3D("BRCA")
  cr.object$plot.type <- "tSNE"
  output <- do.call(tSNE3D, cr.object)

  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  expect3DHtmlOutput(output$plots)
})
