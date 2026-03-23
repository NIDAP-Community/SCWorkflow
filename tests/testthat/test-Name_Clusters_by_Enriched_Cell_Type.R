test_that("Run Name clusters with default parameters - TEC data", {
  input <- getParamsNameClus("TEC")
  output <- do.call(nameClusters, input)

  #Test output values and plot:
  newclus <- output$object@meta.data$clusternames
  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_s4_class(output$object, "Seurat")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawpng(output$plots),
    "TEC_clusters.png"
  )

})

test_that("Run Name clusters with interactive plot", {
  # load data
  input <- getParamsNameClus("TEC")
  input$interactive = TRUE
  output <- do.call(nameClusters, input)

  expect_equal(class(output$plots), c("plotly", "htmlwidget"))
  # expect_snapshot_file(
  #   .drawplot(output$plots),
  #   "TEC_clusters2.png"
  # )
})

test_that("Run Name clusters with ordering celltypes", {
  input <- getParamsNameClus("TEC")
  input$order.celltypes.by = c(
    "B cells",
    "Dendritic cells",
    "Endothelial cells",
    "Macrophages",
    "Monocytes",
    "Epithelial cells",
    "Erythrocytes",
    "Fibroblasts",
    "Hepatocytes",
    "Neurons",
    "T cells"
  )

  output <- do.call(nameClusters, input)

  skip_on_ci()
  expect_snapshot_file(
    .drawpng(output$plots),
    "TEC_clusters_ordered.png"
  )
})

test_that("Run Name clusters with ordering warning missing some celltypes", {
  input <- getParamsNameClus("TEC")
  input$order.celltypes.by = c(
    "B cells",
    "Dendritic cells",
    "Endothelial cells",
    "Macrophages",
    "Monocytes",
    "Epithelial cells",
    #"Erythrocytes","Fibroblasts",
    "Hepatocytes",
    "Neurons",
    "T cells"
  )

  expect_warning(output <- do.call(nameClusters, input),
                 "^Some factors were not included in the list")

  skip_on_ci()
  expect_snapshot_file(
    .drawpng(output$plots),
    "TEC_clusters_missing.png"
  )
})

test_that("Run Name clusters with ordering warning adding unknown celltypes",
          {
            input <- getParamsNameClus("TEC")
            input$order.celltypes.by = c(
              "B cells",
              "Dendritic cells",
              "Endothelial cells",
              "Macrophages",
              "Monocytes",
              "Epithelial cells",
              "celltype1",
              "celltype2",
              "Erythrocytes",
              "Fibroblasts",
              "Hepatocytes",
              "Neurons",
              "T cells"
            )
            expect_warning(output <- do.call(nameClusters, input),
                           "^Some factors are not in data")
            expect_snapshot_file(
              .drawpng(output$plots),
              "TEC_clusters_missing2.png"
            )
          })

test_that("Run Name clusters with default parameters - Chariou", {
  input <- getParamsNameClus("Chariou")
  output <- do.call(nameClusters, input)

  #Test output values and plot:
  newclus <- output$object@meta.data$Clusternames
  expect_equal(sort(unique(newclus)), sort(input$cluster.names))

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_s4_class(output$object, "Seurat")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)

  skip_on_ci()
  expect_snapshot_file(
    .drawpng(output$plots),
    "Chariou_clusters.png"
  )
})

test_that("Run Name clusters with default parameters - PBMC single", {
  input <- getParamsNameClus("pbmc-single")
  output <- do.call(nameClusters, input)

  #Test output values and plot:
  newclus <- output$object@meta.data$Clusternames
  expect_equal(sort(unique(newclus)), sort(input$cluster.names))

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_s4_class(output$object, "Seurat")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)

  skip_on_ci()
  expect_snapshot_file(
    .drawpng(output$plots),
    "NSCLC_single_clusters.png"
  )
})

test_that("Run Name clusters with default parameters - NSCLC multi", {
  input <- getParamsNameClus("nsclc-multi")
  output <- do.call(nameClusters, input)

  #Test output values and plot:
  newclus <- output$object@meta.data$Clusternames
  expect_equal(sort(unique(newclus)), sort(input$cluster.names))

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_s4_class(output$object, "Seurat")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)

  skip_on_ci()
  expect_snapshot_file(
    .drawpng(output$plots),
    "NSCLC_multi_clusters.png"
  )
})
