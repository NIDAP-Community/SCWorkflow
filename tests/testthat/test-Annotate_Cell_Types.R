test_that("Annotate_Cell_Types run with normal parameters - Mouse TEC Data",
          {
            tec.data <- getParamACT("TEC")
            output <- suppressWarnings(do.call(annotateCellTypes, tec.data))
            
            ggsave(
              "output/TEC_annotateCellTypes.p1.png",
              output$p1,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "TEC_annotateCellTypes.p1.png")
            ggsave(
              "output/TEC_annotateCellTypes.p2.png",
              output$p2,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "TEC_annotateCellTypes.p2.png")
            
            expect_type(output, "list")
            expected.elements = c("object", "plots")
            expect_setequal(names(output), expected.elements)
          })

test_that("Annotate_Cell_Types run with reduction type TSNE - Mouse TEC Data",
          {
            tec.data <- getParamACT("TEC")
            tec.data$reduction.type <-  "tsne"
            output <- suppressWarnings(do.call(annotateCellTypes, tec.data))
            
            ggsave(
              "output/TEC_annotateCellTypes.tsne.p1.png",
              output$p1,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "TEC_annotateCellTypes.tsne.p1.png")
            ggsave(
              "output/TEC_annotateCellTypes.tsne.p2.png",
              output$p2,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "TEC_annotateCellTypes.tsne.p2.png")
            
            expect_type(output, "list")
            expected.elements = c("object", "plots")
            expect_setequal(names(output), expected.elements)
          })


test_that("Test Annotate Cell Types with FineTuning - Mouse TEC dataset", {
  tec.data <- getParamACT("TEC")
  tec.data$do.finetuning <- TRUE
  output <- suppressWarnings(do.call(annotateCellTypes, tec.data))
  
  ggsave(
    "output/TEC_annotateCellTypes.fine.p1.png",
    output$p1,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_annotateCellTypes.fine.p1.png")
  ggsave(
    "output/TEC_annotateCellTypes.fine.p2.png",
    output$p2,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "TEC_annotateCellTypes.fine.p2.png")
  
  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})


test_that("Annotate_Cell_Types run with normal parameters - Chariou Data", {
  chariou_data <- getParamACT("Chariou")
  output <- suppressWarnings(do.call(annotateCellTypes, chariou_data))
  
  ggsave(
    "output/Chariou_annotateCellTypes.p1.png",
    output$p1,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "Chariou_annotateCellTypes.p1.png")
  ggsave(
    "output/Chariou_annotateCellTypes.p2.png",
    output$p2,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "Chariou_annotateCellTypes.p2.png")
  
  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
  
})

test_that("Annotate_Cell_Types run with normal parameters - PBMC-single Data",
          {
            pbmc.single.data <- getParamACT("pbmc-single")
            output <- suppressWarnings(do.call(annotateCellTypes, pbmc.single.data))
            
            ggsave(
              "output/PBMCsinlge_annotateCellTypes.p1.png",
              output$p1,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "PBMCsingle_annotateCellTypes.p1.png")
            ggsave(
              "output/PBMCsingle_annotateCellTypes.p2.png",
              output$p2,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "PBMCsingle_annotateCellTypes.p2.png")
            
            expect_type(output, "list")
            expected.elements = c("object", "plots")
            expect_setequal(names(output), expected.elements)
          })

test_that("Annotate_Cell_Types run with normal parameters - NSCLC-multi Data",
          {
            nsclc.multi.data <- getParamACT("nsclc-multi")
            output <- suppressWarnings(do.call(annotateCellTypes, nsclc.multi.data))
            
            ggsave(
              "output/NSCLCmulti_annotateCellTypes.p1.png",
              output$p1,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "NSCLCmulti_annotateCellTypes.p1.png")
            ggsave(
              "output/NSCLCmulti_annotateCellTypes.p2.png",
              output$p2,
              width = 10,
              height = 10
            )
            expect_snapshot_file("output", "NSCLCmulti_annotateCellTypes.p2.png")
            
            expect_type(output, "list")
            expected.elements = c("object", "plots")
            expect_setequal(names(output), expected.elements)
          })

test_that("Annotate_Cell_Types run with normal parameters - BRCA Data", {
  brca.data <- getParamACT("BRCA")
  output <- suppressWarnings(do.call(annotateCellTypes, brca.data))
  
  ggsave(
    "output/BRCA_annotateCellTypes.p1.png",
    output$p1,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "BRCA_annotateCellTypes.p1.png")
  ggsave(
    "output/BRCA_annotateCellTypes.p2.png",
    output$p2,
    width = 10,
    height = 10
  )
  expect_snapshot_file("output", "BRCA_annotateCellTypes.p2.png")
  
  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})
