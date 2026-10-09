test_that("plotMetadata errors when metadata.to.plot is empty", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$metadata.to.plot <- "c()"

  expect_error(
    do.call(plotMetadata, chariou.data),
    "metadata.to.plot.*Available metadata columns:.*SCT_snn_res.2.4"
  )
})

test_that("plotMetadata accepts exact metadata names without renaming object columns", {
  chariou.data <- getParamPM("Chariou")
  original.metadata.columns <- colnames(chariou.data$object@meta.data)
  expect_true("SCT_snn_res.2.4" %in% original.metadata.columns)
  expect_false("SCT_snn_res_2_4" %in% original.metadata.columns)

  chariou.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'
  captured.output <- capture.output(
    output <- do.call(plotMetadata, chariou.data)
  )

  expect_type(output, "list")
  expect_length(output$plots, 1)
  expect_identical(colnames(output$object@meta.data), original.metadata.columns)
  expect_equal(sum(grepl("^Available metadata columns:", captured.output)), 1)
  expect_match(
    paste(captured.output, collapse = "\n"),
    "Available metadata columns:.*SCT_snn_res.2.4"
  )
})

test_that("plotMetadata errors when metadata.to.plot is not an exact column name", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res_2_4")'

  captured.output <- capture.output(
    expect_error(
      do.call(plotMetadata, chariou.data),
      "not found in the object: SCT_snn_res_2_4.*Available metadata columns:.*SCT_snn_res.2.4"
    )
  )
  expect_false(any(grepl("^Available metadata columns:", captured.output)))
})

test_that("plotMetadata accepts dotted and underscore metadata names when both exist", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$object@meta.data$SCT_snn_res_2_4 <-
    chariou.data$object@meta.data[["SCT_snn_res.2.4"]]
  original.metadata.columns <- colnames(chariou.data$object@meta.data)
  expect_true("SCT_snn_res.2.4" %in% original.metadata.columns)
  expect_true("SCT_snn_res_2_4" %in% original.metadata.columns)

  dotted.data <- chariou.data
  dotted.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'
  dotted.captured.output <- capture.output(
    dotted.output <- do.call(plotMetadata, dotted.data)
  )

  underscore.data <- chariou.data
  underscore.data$metadata.to.plot <- 'c("SCT_snn_res_2_4")'
  underscore.captured.output <- capture.output(
    underscore.output <- do.call(plotMetadata, underscore.data)
  )

  expect_length(dotted.output$plots, 1)
  expect_identical(colnames(dotted.output$object@meta.data), original.metadata.columns)
  expect_equal(sum(grepl("^Available metadata columns:", dotted.captured.output)), 1)
  expect_length(underscore.output$plots, 1)
  expect_identical(colnames(underscore.output$object@meta.data), original.metadata.columns)
  expect_equal(sum(grepl("^Available metadata columns:", underscore.captured.output)), 1)
})

test_that("plotMetadata accepts orig_ident without renaming metadata columns", {
  chariou.data <- getParamPM("Chariou")
  colnames(chariou.data$object@meta.data)[colnames(chariou.data$object@meta.data) == "orig.ident"] <-
    "orig_ident"
  chariou.data$samples.to.include <- "c()"
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'
  original.metadata.columns <- colnames(chariou.data$object@meta.data)
  expect_true("orig_ident" %in% original.metadata.columns)
  expect_false("orig.ident" %in% original.metadata.columns)

  captured.output <- capture.output(
    output <- do.call(plotMetadata, chariou.data)
  )

  expect_length(output$plots, 1)
  expect_identical(colnames(output$object@meta.data), original.metadata.columns)
  expect_equal(sum(grepl("^Available metadata columns:", captured.output)), 1)
})

test_that("Test Plot Metadata using TEC (Mouse) dataset", {
  tec.data <- getParamPM("TEC")
  output <- do.call(plotMetadata,tec.data)

  ggsave("output/TEC_plotmet.png",output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet1.png")
  ggsave("output/TEC_plotmet.png",output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet2.png")
  ggsave("output/TEC_plotmet.png",output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet3.png")
  ggsave("output/TEC_plotmet.png",output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet4.png")
  ggsave("output/TEC_plotmet.png",output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet5.png")
  ggsave("output/TEC_plotmet.png",output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet6.png")
  
    
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
  

})


test_that("Test Plot Metadata using Chariou (Mouse) dataset", {
  
  chariou.data <- getParamPM("Chariou")
  output <- do.call(plotMetadata,chariou.data)

  ggsave("output/Chariou_plotmet.png",output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output","Chariou_plotmet1.png")
  ggsave("output/Chariou_plotmet.png",output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output","Chariou_plotmet2.png")
  ggsave("output/Chariou_plotmet.png",output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output","Chariou_plotmet3.png")
  
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
  
})



test_that("Test Plot Metadata using BRCA (Human) dataset", {

  brca.data <- getParamPM("BRCA")
  output <- do.call(plotMetadata, brca.data)

  ggsave("output/BRCA_plotmet.png",output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output","BRCA_plotmet1.png")
  ggsave("output/BRCA_plotmet.png",output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output","BRCA_plotmet2.png")
  ggsave("output/BRCA_plotmet.png",output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output","BRCA_plotmet3.png")
  ggsave("output/BRCA_plotmet.png",output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output","BRCA_plotmet4.png")
  ggsave("output/BRCA_plotmet.png",output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output","BRCA_plotmet5.png")
  ggsave("output/BRCA_plotmet.png",output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output","BRCA_plotmet6.png")
  
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)

})




test_that("Test Plot Metadata using NSCLCmulti (Human) dataset", {
  
  nsclc.multi.data <- getParamPM("nsclc-multi")
  output <- do.call(plotMetadata,nsclc.multi.data)

  ggsave("output/NSCLCmulti_plotmet.png",output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output","NSCLCmulti_plotmet1.png")
  ggsave("output/NSCLCmulti_plotmet.png",output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output","NSCLCmulti_plotmet2.png")
  ggsave("output/NSCLCmulti_plotmet.png",output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output","NSCLCmulti_plotmet3.png")
  ggsave("output/NSCLCmulti_plotmet.png",output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output","NSCLCmulti_plotmet4.png")
  ggsave("output/NSCLCmulti_plotmet.png",output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output","NSCLCmulti_plotmet5.png")
  ggsave("output/NSCLCmulti_plotmet.png",output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output","NSCLCmulti_plotmet6.png")
  
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)

})



test_that("Test Plot Metadata using PBMCsingle (Human) dataset", {


  pbmc.single.data <- getParamPM("pbmc-single")
  output <- do.call(plotMetadata,pbmc.single.data)

  ggsave("output/PBMCsingle_plotmet.png",output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output","PBMCsingle_plotmet1.png")
  ggsave("output/PBMCsingle_plotmet.png",output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output","PBMCsingle_plotmet2.png")
  ggsave("output/PBMCsingle_plotmet.png",output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output","PBMCsingle_plotmet3.png")
  ggsave("output/PBMCsingle_plotmet.png",output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output","PBMCsingle_plotmet4.png")
  ggsave("output/PBMCsingle_plotmet.png",output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output","PBMCsingle_plotmet5.png")
  ggsave("output/PBMCsingle_plotmet.png",output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output","PBMCsingle_plotmet6.png")
  
  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
  
  
})




test_that("Test Plot Metadata using TEC (Mouse) dataset; UMAP", {


  tec.data <- getParamPM("TEC")
  tec.data$reduction.type <- "umap"

  output <- do.call(plotMetadata,tec.data)

  ggsave("output/TEC_plotmet.png",output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet1.umap.png")
  ggsave("output/TEC_plotmet.png",output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet2.umap.png")
  ggsave("output/TEC_plotmet.png",output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet3.umap.png")
  ggsave("output/TEC_plotmet.png",output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet4.umap.png")
  ggsave("output/TEC_plotmet.png",output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet5.umap.png")
  ggsave("output/TEC_plotmet.png",output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output","TEC_plotmet6.umap.png")

  expect_type(output,"list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
  

})
