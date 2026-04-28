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
