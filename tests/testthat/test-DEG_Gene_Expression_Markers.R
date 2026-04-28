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
  output <- do.call(degGeneExpressionMarkers, chariou.data)
  
  expect_type(output, "list")
  expected.elements = c("data")
  expect_setequal(names(output), expected.elements)
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
