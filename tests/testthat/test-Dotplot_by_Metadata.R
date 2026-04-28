test_that("dotplot run with normal parameters - TEC Data", {
  tec.data <- getParamDP("TEC")
  output <- do.call(dotPlotMet, tec.data)

  expect_type(output, "list")
  expected.elements <- c("plots","data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawDotPng(output$plot),
    "TEC_dotplot.png"
  )
})

test_that("dotplot run with reversed axes - TEC Data", {
  tec.data <- getParamDP("TEC")
  tec.data$plot.reverse <- TRUE
  output <- do.call(dotPlotMet, tec.data)
  
  expect_type(output, "list")
  expected.elements <- c("plots","data")
  expect_setequal(names(output), expected.elements)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawDotPng(output$plot),
    "TEC_dotplot_reversed.png"
  )
})

test_that("dotplot run with different color - TEC Data", {
  tec.data <- getParamDP("TEC")
  tec.data$dot.color <- "red"
  output <- do.call(dotPlotMet, tec.data)

  expect_type(output, "list")
  expected.elements <- c("plots","data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawDotPng(output$plot),
    "TEC_dotplot_red.png"
  )
})

test_that("dotplot run with warning for fewer categories - TEC Data", {
  tec.data <- getParamDP("TEC")
  tec.data$cells <- c(0, 1, 2, 3, 4, 5, 6)
  expect_warning(do.call(dotPlotMet, tec.data),
                 "^There are")
})

test_that("dotplot run with error for non-matched category - TEC Data", {
  tec.data <- getParamDP("TEC")
  tec.data <- tec.data
  tec.data$metadata <- "orig.ident"
  expect_error(
    do.call(dotPlotMet, tec.data),
    "^At least 2 metadata categories you wish to plot should"
  )
})

test_that("dotplot run with warning for duplicate genes", {
  tec.data <- getParamDP("TEC")
  tec.data <- tec.data
  tec.data$markers <- c("Map7", tec.data$markers)
  expect_warning(
    do.call(dotPlotMet, tec.data),
    fixed = TRUE,
    "The following duplicate genes were removed: Map7"
  )
})

test_that("dotplot run with warning for missing genes", {
  tec.data <- getParamDP("TEC")
  tec.data <- tec.data
  tec.data$markers <- c("Adgre1", "Ccr2", tec.data$markers)
  expect_warning(
    do.call(dotPlotMet, tec.data),
    "^There are "
  )

  #Full message:
  #There are 2 gene(s) absent from dataset:'Adgre1', 'Ccr2'.
  #Possible reasons are that gene is not official gene symbol or
  #gene is not highly expressed and has been filtered.
})

test_that("dotplot run with error for missing all genes", {
  tec.data <- getParamDP("TEC")
  tec.data <- tec.data
  tec.data$markers <- c("APOE", "ARG1", "CD38", "CD3D", "CD3E", "CD3G")
  expect_error(do.call(dotPlotMet, tec.data),
               "No genes listed are found in dataset.")

})

test_that("dotplot produced and contingency table returned - Chariou Data", {
  chariou.data <- getParamDP("Chariou")
  output <- do.call(dotPlotMet, chariou.data)

  expect_type(output, "list")
  expected.elements <- c("plots","data")
  expect_setequal(names(output), expected.elements)
  
  skip_on_ci()
  expect_snapshot_file(
      .drawDotPng(output$plot),
      "Chariou_dotplot.png"
  )
})

test_that("dotplot produced and contingency table returned - PBMC-single Data",
    {
      pbmc.single.data <- getParamDP("pbmc-single")
        output <- do.call(dotPlotMet, pbmc.single.data)

        expect_type(output, "list")
        expected.elements <- c("plots","data")
        expect_setequal(names(output), expected.elements)

        skip_on_ci()
        expect_snapshot_file(
          .drawDotPng(output$plot),
          "pbmc-single_dotplot.png"
        )
})

test_that("dotplot produced and contingency table returned - NSCLC-multi Data",
    {
      nsclc.multi.data <- getParamDP("nsclc-multi")
      expect_warning(output <- do.call(dotPlotMet, nsclc.multi.data),
                     "^There are")

      expect_type(output, "list")
      expected.elements <- c("plots","data")
      expect_setequal(names(output), expected.elements)

      skip_on_ci()
      expect_snapshot_file(
        .drawDotPng(output$plot),
        "nsclc-multi_dotplot.png"
      )
})

test_that("dotplot produced and contingency table returned - BRCA Data", {
  brca.data <- getParamDP("BRCA")
  expect_warning(output <- do.call(dotPlotMet, brca.data),
                 "^There are")

  expect_type(output, "list")
  expected.elements <- c("plots","data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(
    .drawDotPng(output$plot),
    "brca_dotplot.png"
  )
})
