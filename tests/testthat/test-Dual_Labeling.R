test_that("Test Dual labeling TEC Data", {
  CRObject <- getParamDL("TEC")
  output <- do.call(dualLabeling, CRObject)
  
  expected.elements = c("object","data" ,"plots")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plots$tsne, 6)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plots$tsne),
        "TEC_duallabel.png"
    )
  expect_snapshot_file(
    .drawdualplot(output$plots$densityHM),
    "TEC_duallabelDensity.png"
  )
  expect_snapshot_file(
    .drawdualtable(output$data$plot_table),
    "TEC_dualtable.png"
  )
})

test_that("Dual labeling with missing gene", {
  CRObject <- getParamDL("TEC")
  CRObject$marker.1 <- "Cd8"
  expect_error(do.call(dualLabeling, CRObject),
               "is not found in dataset$")
})

test_that("Dual labeling with umap", {
  CRObject <- getParamDL("TEC")
  CRObject$data.reduction = "umap"
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object","data" ,"plots")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plots$umap, 6)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plots$umap),
             "umap_duallabel.png"
    )
})


test_that("Test Dual labeling Chariou Data", {
  CRObject <- getParamDL("Chariou")
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object","data" ,"plots")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plots$tsne, 6)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plots$tsne),
       "chariou_duallabel.png"
    )
  expect_snapshot_file(
    .drawdualplot(output$plots$densityHM),
    "chariou_duallabelDensity.png"
  )
  expect_snapshot_file(
    .drawdualtable(output$data$plot_table),
    "chariou_dualtable.png"
  )
})

test_that("Test Dual labeling NSCLC_Single Data", {
  CRObject <- getParamDL("NSCLC_Single")
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object","data" ,"plots")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plots$tsne, 6)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plots$tsne),
             "NSCLC_Single_duallabel.png"
    )
  expect_snapshot_file(
    .drawdualplot(output$plots$densityHM),
    "NSCLC_Single_duallabelDensity.png"
  )
  expect_snapshot_file(
    .drawdualtable(output$data$plot_table),
    "NSCLC_Single_dualtable.png"
  )
})

test_that("Test Dual labeling NSCLC-multi Data", {
  CRObject <- getParamDL("nsclc-multi")
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object","data" ,"plots")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plots$tsne, 6)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plots$tsne),
             "nsclc-multi_duallabel.png"
    )
  expect_snapshot_file(
    .drawdualplot(output$plots$densityHM),
    "nsclc-multi_duallabelDensity.png"
  )
  expect_snapshot_file(
    .drawdualtable(output$data$plot_table),
    "nsclc-multi_dualtable.png"
  )
})

test_that("Test Dual labeling BRCA Data", {
  CRObject <- getParamDL("BRCA")
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object","data" ,"plots")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plots$tsne, 6)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plots$tsne),
             "brca_duallabel.png"
    )
  expect_snapshot_file(
    .drawdualplot(output$plots$densityHM),
    "brca_duallabelDensity.png"
  )
  expect_snapshot_file(
    .drawdualtable(output$data$plot_table),
    "brca_dualtable.png"
  )
})
