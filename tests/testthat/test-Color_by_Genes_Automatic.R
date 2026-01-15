# Data Testing
test_that("Color by Genes Automatic works for TEC data", {

  tec.data <- getCbgAutoParam("TEC")

  cbg.demo <- colorByMarkerTable(object = tec.data$object,
                           samples.subset = tec.data$samples.subset,
                           samples.to.display = tec.data$samples.to.display,
                           marker.table = tec.data$marker.table,
                           cells.of.interest = tec.data$cells.of.interest)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCbG(cbg.demo),
    "tec_cbg.png"
  )

  expect_true(is.list(cbg.demo))
  expect_setequal(names(cbg.demo), c("overall", "celltype", "manual_entry"))
  expect_true(is.list(cbg.demo$overall))
  expect_true(is.list(cbg.demo$celltype))
  expect_gt(length(cbg.demo$overall), 0)
  expect_gt(length(cbg.demo$celltype), 0)
  expect_true(grid::is.grob(cbg.demo$overall[[1]]))
  expect_true(grid::is.grob(cbg.demo$celltype[[1]]))
  expect_true(is.null(cbg.demo$manual_entry) || grid::is.grob(cbg.demo$manual_entry))
})

test_that("Color by Genes Automatic works for Chariou data", {

  chariou.data <- getCbgAutoParam("Chariou")

  cbg.demo <- colorByMarkerTable(object = chariou.data$object,
                           samples.subset = chariou.data$samples.subset,
                           samples.to.display = chariou.data$samples.to.display,
                           marker.table = chariou.data$marker.table,
                           cells.of.interest = chariou.data$cells.of.interest)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCbG(cbg.demo),
    "chariou_cbg.png"
  )

  expect_true(is.list(cbg.demo))
  expect_setequal(names(cbg.demo), c("overall", "celltype", "manual_entry"))
  expect_true(is.list(cbg.demo$overall))
  expect_true(is.list(cbg.demo$celltype))
  expect_gt(length(cbg.demo$overall), 0)
  expect_gt(length(cbg.demo$celltype), 0)
  expect_true(grid::is.grob(cbg.demo$overall[[1]]))
  expect_true(grid::is.grob(cbg.demo$celltype[[1]]))
  expect_true(is.null(cbg.demo$manual_entry) || grid::is.grob(cbg.demo$manual_entry))

})

test_that("Color by Genes Automatic works for pbmc.single data", {

  pbmc.single <- getCbgAutoParam("pbmc.single")

  cbg.demo <- colorByMarkerTable(object = pbmc.single$object,
                           samples.subset = pbmc.single$samples.subset,
                           samples.to.display = pbmc.single$samples.to.display,
                           marker.table = pbmc.single$marker.table,
                           cells.of.interest = pbmc.single$cells.of.interest)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCbG(cbg.demo),
    "pbmc_single_cbg.png"
  )

  expect_true(is.list(cbg.demo))
  expect_setequal(names(cbg.demo), c("overall", "celltype", "manual_entry"))
  expect_true(is.list(cbg.demo$overall))
  expect_true(is.list(cbg.demo$celltype))
  expect_gt(length(cbg.demo$overall), 0)
  expect_gt(length(cbg.demo$celltype), 0)
  expect_true(grid::is.grob(cbg.demo$overall[[1]]))
  expect_true(grid::is.grob(cbg.demo$celltype[[1]]))
  expect_true(is.null(cbg.demo$manual_entry) || grid::is.grob(cbg.demo$manual_entry))

})

test_that("Color by Genes Automatic works for nsclc_multi data", {

  nsclc_multi <- getCbgAutoParam("nsclc.multi")

  cbg.demo <- colorByMarkerTable(object = nsclc_multi$object,
                           samples.subset = nsclc_multi$samples.subset,
                           samples.to.display = nsclc_multi$samples.to.display,
                           marker.table = nsclc_multi$marker.table,
                           cells.of.interest = nsclc_multi$cells.of.interest)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCbG(cbg.demo),
    "nsclc_multi_cbg.png"
  )

  expect_true(is.list(cbg.demo))
  expect_setequal(names(cbg.demo), c("overall", "celltype", "manual_entry"))
  expect_true(is.list(cbg.demo$overall))
  expect_true(is.list(cbg.demo$celltype))
  expect_gt(length(cbg.demo$overall), 0)
  expect_gt(length(cbg.demo$celltype), 0)
  expect_true(grid::is.grob(cbg.demo$overall[[1]]))
  expect_true(grid::is.grob(cbg.demo$celltype[[1]]))
  expect_true(is.null(cbg.demo$manual_entry) || grid::is.grob(cbg.demo$manual_entry))

})

test_that("Color by Genes Automatic works for BRCA data", {

  BRCA_data <- getCbgAutoParam("BRCA")

  cbg.demo <- colorByMarkerTable(object = BRCA_data$object,
                           samples.subset = BRCA_data$samples.subset,
                           samples.to.display = BRCA_data$samples.to.display,
                           marker.table = BRCA_data$marker.table,
                           cells.of.interest = BRCA_data$cells.of.interest)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCbG(cbg.demo),
    "brca_cbg.png"
  )

  expect_true(is.list(cbg.demo))
  expect_setequal(names(cbg.demo), c("overall", "celltype", "manual_entry"))
  expect_true(is.list(cbg.demo$overall))
  expect_true(is.list(cbg.demo$celltype))
  expect_gt(length(cbg.demo$overall), 0)
  expect_gt(length(cbg.demo$celltype), 0)
  expect_true(grid::is.grob(cbg.demo$overall[[1]]))
  expect_true(grid::is.grob(cbg.demo$celltype[[1]]))
  expect_true(is.null(cbg.demo$manual_entry) || grid::is.grob(cbg.demo$manual_entry))

})


# Error Checking
test_that("Color by Genes Automatic stops when user inputs an assay not found
          in seurat ", {

  tec.data <- getCbgAutoParam("TEC")

  expect_error(colorByMarkerTable(object = tec.data$object,
                            samples.subset = tec.data$samples.subset,
                            samples.to.display = tec.data$samples.to.display,
                            marker.table = tec.data$marker.table,
                            cells.of.interest = tec.data$cells.of.interest,
                            assay = "wrong_assay"))

})

test_that("Color by Genes Automatic stops when user inputs an reduction type
          not found in seurat ", {

  tec.data <- getCbgAutoParam("TEC")

  expect_error(colorByMarkerTable(object = tec.data$object,
                            samples.subset = tec.data$samples.subset,
                            samples.to.display = tec.data$samples.to.display,
                            marker.table = tec.data$marker.table,
                            cells.of.interest = tec.data$cells.of.interest,
                            reduction.type = "wrong_reduction"),
                            NA)

})
