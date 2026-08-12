test_that(
  "Test Filter Seurat Object by Metadata using Downsampled TEC (Mouse) data
   with normal parameters",
  {
    tec.data <- getParamFSOBM("TEC")
    output <- do.call(filterSeuratObjectByMetadata, tec.data)
    
    ggsave("output/TEC_fsobm.plot1.png",
           output$plot1,
           width = 10,
           height = 10)
    expect_snapshot_file("output", "TEC_fsobm.plot1.png")
    ggsave("output/TEC_fsobm.plot2.png",
           output$plot2,
           width = 10,
           height = 10)
    expect_snapshot_file("output", "TEC_fsobm.plot2.png")
    
    expect_type(output, "list")
    expected.elements = c("object", "plots")
    expect_setequal(names(output), expected.elements)
    
  }
)


test_that("Test Filter Seurat Object by Metadata using Chariou (Mouse) dataset",
          {
            chariou.data <- getParamFSOBM("Chariou")
            output <-
              do.call(filterSeuratObjectByMetadata, chariou.data)
            
            ggsave(
              "output/Chariou_fsobm.plot1.png",
              output$plot1,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "Chariou_fsobm.plot1.png")
            ggsave(
              "output/Chariou_fsobm.plot2.png",
              output$plot2,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "Chariou_fsobm.plot2.png")
            
            expect_type(output, "list")
            expected.elements = c("object", "plots")
            expect_setequal(names(output), expected.elements)
            
          })

test_that("Blank sample name column uses all samples", {
  expect_identical(deparse(formals(filterSeuratObjectByMetadata)$samples.to.include), "c(\"\")")
  expect_identical(formals(filterSeuratObjectByMetadata)$sample.name, "")

  chariou.data <- getParamFSOBM("Chariou")
  chariou.data$sample.name <- ""
  chariou.data$samples.to.include <- ""
  output <- do.call(filterSeuratObjectByMetadata, chariou.data)

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
  expect_true(ncol(output$object) > 0)
})

test_that("Dotted metadata category names are matched after normalization", {
  chariou.data <- getParamFSOBM("Chariou")
  chariou.data$object@meta.data[["SCT_snn_res.0.4"]] <- chariou.data$object@meta.data$seurat_clusters
  chariou.data$sample.name <- ""
  chariou.data$samples.to.include <- c("")
  chariou.data$category.to.filter <- "SCT_snn_res.0.4"
  chariou.data$values.to.filter <- c("0", "1")
  output <- do.call(filterSeuratObjectByMetadata, chariou.data)

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
  expect_true(ncol(output$object) > 0)
})



test_that("Test Filter Seurat Object by Metadata using BRCA (Human) dataset",
          {
            brca.data <- getParamFSOBM("BRCA")
            output <-
              do.call(filterSeuratObjectByMetadata, brca.data)
            
            ggsave(
              "output/BRCA_fsobm.plot1.png",
              output$plot1,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "BRCA_fsobm.plot1.png")
            ggsave(
              "output/BRCA_fsobm.plot2.png",
              output$plot2,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "BRCA_fsobm.plot2.png")
            
            expect_type(output, "list")
            expected.elements = c("object", "plots")
            expect_setequal(names(output), expected.elements)
            
          })


test_that("Test Filter Seurat Object by Metadata using NSCLCmulti (Human) data",
          {
            nsclc.multi.data <- getParamFSOBM("nsclc-multi")
            output <-
              do.call(filterSeuratObjectByMetadata, nsclc.multi.data)
            
            ggsave(
              "output/NSCLCmulti_fsobm.plot1.png",
              output$plot1,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "NSCLCmulti_fsobm.plot1.png")
            ggsave(
              "output/NSCLCmulti_fsobm.plot2.png",
              output$plot2,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "NSCLCmulti_fsobm.plot2.png")
            
            expect_type(output, "list")
            expected.elements = c("object", "plots")
            expect_setequal(names(output), expected.elements)
            
          })




test_that("Test Filter Seurat Object by Metadata using PBMCsingle (Human) data",
          {
            pbmc.single.data <- getParamFSOBM("pbmc-single")
            output <-
              do.call(filterSeuratObjectByMetadata, pbmc.single.data)
            
            ggsave(
              "output/PBMCsingle_fsobm.plot1.png",
              output$plot1,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "PBMCsingle_fsobm.plot1.png")
            ggsave(
              "output/PBMCsingle_fsobm.plot2.png",
              output$plot2,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "PBMCsingle_fsobm.plot2.png")
            
            expect_type(output, "list")
            expected.elements = c("object", "plots")
            expect_setequal(names(output), expected.elements)
            
          })
