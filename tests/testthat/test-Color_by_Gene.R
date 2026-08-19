test_that("Color by Gene using TEC (Mouse) dataset with normal parameters", {
  tec.data <- getParamCBG("TEC")
  output <- do.call(colorByGene, tec.data)
  
  ggsave(
    "output/TEC_colbygeneplot_Gapdh.png",
    output$plots[[1]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_colbygeneplot_Gapdh.png")
  
  ggsave(
    "output/TEC_colbygeneplot_Il2ra.png",
    output$plots[[2]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_colbygeneplot_Il2ra.png")
  
  expect_type(output, "list")
  expected.elements <- c("data", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Color by Gene using Chariou (Mouse) dataset", {
  chariou.data <- getParamCBG("Chariou")
  output <- do.call(colorByGene, chariou.data)
  
  ggsave(
    "output/Chariou_colbygeneplot_Gapdh.png",
    output$plots[[1]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "Chariou_colbygeneplot_Gapdh.png")
  
  expect_type(output, "list")
  expected.elements <- c("data", "plots")
  expect_setequal(names(output), expected.elements)
})


test_that("Color by Gene using BRCA (Human) dataset", {
  brca.data <- getParamCBG("BRCA")
  output <- do.call(colorByGene, brca.data)
  
  ggsave(
    "output/BRCA_colbygeneplot_PARP8.png",
    output$plots[[1]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "BRCA_colbygeneplot_PARP8.png")
  
  expect_type(output, "list")
  expected.elements <- c("data", "plots")
  expect_setequal(names(output), expected.elements)
})


test_that("Test Annotate Cell Types using NSCLCmulti (Human) dataset", {
  nsclc.multi.data <- getParamCBG("nsclc-multi")
  output <- do.call(colorByGene, nsclc.multi.data)
  
  ggsave(
    "output/NSCLCmulti_colbygeneplot_AMER2.png",
    output$plots[[1]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "NSCLCmulti_colbygeneplot_AMER2.png")

  ggsave(
    "output/NSCLCmulti_colbygeneplot_CD9.png",
    output$plots[[2]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "NSCLCmulti_colbygeneplot_CD9.png")
  
  ggsave(
    "output/NSCLCmulti_colbygeneplot_CEBPD.png",
    output$plots[[3]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "NSCLCmulti_colbygeneplot_CEBPD.png")
  
  expect_type(output, "list")
  expected.elements <- c("data", "plots")
  expect_setequal(names(output), expected.elements)
})



test_that("Test Annotate Cell Types using PBMCsingle (Human) dataset", {
  pbmc.single.data <- getParamCBG("pbmc-single")
  output <- do.call(colorByGene, pbmc.single.data)
  
  ggsave(
    "output/PBMCsingle_colbygeneplot_MYC.png",
    output$plots[[1]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "PBMCsingle_colbygeneplot_MYC.png")
  
  expect_type(output, "list")
  expected.elements <- c("data", "plots")
  expect_setequal(names(output), expected.elements)
})


test_that("Color by Gene using TEC (Mouse) dataset, TSNE", {
  tec.data <- getParamCBG("TEC")
  tec.data$reduction.type = "tsne"
  
  output <- do.call(colorByGene, tec.data)
  
  ggsave(
    "output/TEC_colbygeneplot_Gapdh.tsne.png",
    output$plots[[1]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_colbygeneplot_Gapdh.tsne.png")
  
  ggsave(
    "output/TEC_colbygeneplot_Il2ra.tsne.png",
    output$plots[[2]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_colbygeneplot_Il2ra.tsne.png")
  
  expect_type(output, "list")
  expected.elements <- c("data", "plots")
  expect_setequal(names(output), expected.elements)
})


test_that("Color by Gene using TEC (Mouse) dataset, Blue color", {
  tec.data <- getParamCBG("TEC")
  tec.data$color <- "blue"
  
  output <- do.call(colorByGene, tec.data)
  
  ggsave(
    "output/TEC_colbygeneplot_Gapdh.tsne.png",
    output$plots[[1]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_colbygeneplot_Gapdh.blue.png")
  
  ggsave(
    "output/TEC_colbygeneplot_Il2ra.tsne.png",
    output$plots[[2]],
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_colbygeneplot_Il2ra.blue.png")
  
  expect_type(output, "list")
  expected.elements <- c("data", "plots")
  expect_setequal(names(output), expected.elements)
})
