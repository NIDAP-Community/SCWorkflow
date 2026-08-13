test_that("Test Plot Metadata using TEC (Mouse) dataset with normal parameters",
          {
            tec.data <- getParamDGEM("TEC")
            output <- do.call(degGeneExpressionMarkers, tec.data)
            
            expect_type(output, "list")
            expected.elements = c("data")
            expect_setequal(names(output), expected.elements)
          })


test_that("Test DEG Gene Expression Markers using negbinom (TEC Mouse dataset)",
          {
            tec.data <- getParamDGEM("TEC")
            tec.data$test.to.use = "negbinom"
            
            output <- do.call(degGeneExpressionMarkers, tec.data)
            
            expect_type(output, "list")
            expected.elements = c("data")
            expect_setequal(names(output), expected.elements)
          })


test_that("Test DEG Gene Expression Markers using Chariou (Mouse) dataset", {
  chariou.data <- getParamDGEM("Chariou")
  logs <- capture.output(output <- do.call(degGeneExpressionMarkers, chariou.data))
  
  expect_type(output, "list")
  expected.elements = c("data")
  expect_setequal(names(output), expected.elements)
  expect_true(any(grepl("Possible values for SCT_snn_res_2_8:", logs)))
  expect_true(any(grepl("  - 0", logs, fixed = TRUE)))
})


test_that("Blank samples uses all samples", {
  chariou.data <- getParamDGEM("Chariou")
  chariou.data$samples <- ""
  logs <- capture.output(output <- do.call(degGeneExpressionMarkers, chariou.data))
  
  expect_type(output, "list")
  expected.elements = c("data")
  expect_setequal(names(output), expected.elements)
  expect_gt(nrow(output$data$DEG_Table), 0)
  expect_true(any(grepl("Possible sample names:", logs)))
  expect_true(any(grepl("CD8dep", logs)))
})


test_that("Omitted samples uses all samples", {
  chariou.data <- getParamDGEM("Chariou")
  chariou.data$samples <- NULL
  output <- do.call(degGeneExpressionMarkers, chariou.data)
  
  expect_type(output, "list")
  expected.elements = c("data")
  expect_setequal(names(output), expected.elements)
  expect_gt(nrow(output$data$DEG_Table), 0)
})


test_that("Omitted parameter prints discrete metadata columns and defaults", {
  chariou.data <- getParamDGEM("Chariou")
  chariou.data$parameter.to.test <- NULL
  logs <- capture.output(output <- do.call(degGeneExpressionMarkers, chariou.data))
  
  expect_type(output, "list")
  expected.elements = c("data")
  expect_setequal(names(output), expected.elements)
  expect_gt(nrow(output$data$DEG_Table), 0)
  expect_true(any(grepl("Possible parameter.to.test columns with categorical/discrete values:", logs)))
  expect_true(any(grepl("SCT_snn_res_2_4", logs)))
  expect_true(any(grepl("No parameter selected, defaulting to", logs)))
  expect_true(any(grepl("Possible values for", logs)))
})


test_that("Test DEG Gene Expression Markers using BRCA (Human) dataset", {
  brca.data <- getParamDGEM("BRCA")
  output <- do.call(degGeneExpressionMarkers, brca.data)
  
  expect_type(output, "list")
  expected.elements = c("data")
  expect_setequal(names(output), expected.elements)
})



test_that("Test DEG Gene Expression Markers using NSCLCmulti (Human) dataset",
          {
            nsclc.multi.data <- getParamDGEM("nsclc-multi")
            output <- do.call(degGeneExpressionMarkers, nsclc.multi.data)
            
            expect_type(output, "list")
            expected.elements = c("data")
            expect_setequal(names(output), expected.elements)
          })



test_that("Test DEG Gene Expression Markers using PBMCsingle (Human) dataset",
          {
            pbmc.single.data <- getParamDGEM("pbmc-single")
            output <- do.call(degGeneExpressionMarkers, pbmc.single.data)
            
            expect_type(output, "list")
            expected.elements = c("data")
            expect_setequal(names(output), expected.elements)
          })
