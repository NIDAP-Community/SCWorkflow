test_that("Produce heatmap and return plot and filtered dataframe: TEC data",
{
  cr.object <- getParamHM("TEC")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)
  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                                 "TEC_heatmap.png")
})

test_that("Heatmap scaled vs unscaled", {
  cr.object <- getParamHM("TEC")
  output <- do.call(heatmapSC, cr.object)
  cr.object$scale.data <- FALSE
  output2 <- do.call(heatmapSC, cr.object)

  #compare scaled (a) vs. nonscaled data (b) to be different
  a <- rowMeans(as.data.frame.matrix(output$data)[, -1])
  b <- rowMeans(as.data.frame.matrix(output2$data)[, -1])

  expect_false(isTRUE(all.equal(a, b)))

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output2$plots),
                       "TEC_heatmap_unscaled.png")

})

test_that("Heatmap run with bad gene name", {
  cr.object <- getParamHM("TEC")

  #Spike in with bad gene name
  cr.object$transcripts <- c("badgene", cr.object$transcripts)

  expect_warning(do.call(heatmapSC, cr.object), "^There are")

  #Actual warning:
  #There are 1 gene(s) absent from dataset:'badgene'.
  #Possible reasons are that gene is not official gene symbol
  #or gene is not highly expressed and has been filtered.

})

test_that("Heatmap run with duplicate gene names", {
  cr.object <- getParamHM("TEC")

  #Spike in with bad gene name
  cr.object$transcripts <- c("Map7", cr.object$transcripts)

  expect_warning(do.call(heatmapSC, cr.object),
                 "The following duplicate genes were removed: Map7")

})

test_that("Heatmap run with error for missing all genes", {
  cr.object <- getParamHM("TEC")
  cr.object$transcripts <-
    c("APOE", "ARG1", "CD38", "CD3D", "CD3E", "CD3G")
  expect_error(do.call(heatmapSC, cr.object),
               "No genes listed are found in dataset.")

})

test_that("Produce heatmap - Chariou data", {
  cr.object <- getParamHM("Chariou")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "Chariou_heatmap.png")
})

test_that("Chariou with no additional protein/transcript annotations", {
  cr.object <- getParamHM("Chariou2")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "Chariou_heatmap2.png")
})

test_that("Produce heatmap - Chariou with no transcripts/proteins", {
  cr.object <- getParamHM("Chariou3")

  expect_error(do.call(heatmapSC, cr.object),
               "At least 1 transcript and/or protein is needed for plotting")
})

test_that("Produce heatmap - PBMC single data", {
  cr.object <- getParamHM("pbmc-single")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "pbmc-single_heatmap.png")
})

test_that("Produce heatmap with protein using NSCLC multi data", {
  cr.object <- getParamHM("nsclc-multi")
  cr.object$legend.font = 3
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "nsclc-multi_heatmap.png")
})

test_that("Produce heatmap with filtered dataframe - BRCA data", {
  cr.object <- getParamHM("BRCA")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "BRCA_heatmap.png")
})