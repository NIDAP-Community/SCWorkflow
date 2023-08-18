# Dataset Testing
test_that("modScore returns metadata with scores and cell calls for tec", {

  tec <- getModuleScoreParam("tec")

  suppressWarnings(modscore.demo <- do.call(modScore, tec))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "tec_MS_1.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "tec_MS_2.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "tec_MS_3.png"
  )

  saveRDS(modscore.demo$harm.object,
          file = test_path("output/modScore/TEC_so.rds"))

  expect_equal_to_reference(modscore.demo$harm.object,
                            file = test_path("output/modScore/TEC_so.rds"))

  expected_elements <- c("Likely_CellType",tec$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
                          failure_message = "modscore results not found")

})

test_that("modScore returns metadata with scores and cell calls for chariou", {

  chariou <- getModuleScoreParam("chariou")

  suppressWarnings(modscore.demo <- do.call(modScore, chariou))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "chariou_MS_1.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "chariou_MS_2.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "chariou_MS_3.png"
  )


  saveRDS(modscore.demo$harm.object,
          file = test_path("output/modScore/chariou_so.rds"))

  expect_equal_to_reference(modscore.demo$harm.object,
                            file = test_path("output/modScore/chariou_so.rds"))

  expected_elements <- c("Likely_CellType",chariou$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
         failure_message = "modscore results not found")

})

test_that("modScore returns metadata with scores and cell calls for
          pbmc.single", {

  pbmc.single <- getModuleScoreParam("pbmc.single")

  suppressWarnings(modscore.demo <- do.call(modScore, pbmc.single))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "pbmc_single_MS_1.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "pbmc_single_MS_2.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "pbmc_single_MS_3.png"
  )

  saveRDS(modscore.demo$harm.object,
          file = test_path("output/modScore/pbmc_single_so.rds"))

  expect_equal_to_reference(modscore.demo$harm.object, file =
                              test_path("output/modScore/pbmc_single_so.rds"))

  expected_elements <- c("Likely_CellType",pbmc.single$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
         failure_message = "modscore results not found")

})

test_that("modScore returns metadata with scores and cell calls for
          nsclc.multi", {

  nsclc.multi <- getModuleScoreParam("nsclc.multi")

  suppressWarnings(modscore.demo <- do.call(modScore, nsclc.multi))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "nsclc_multi_MS_1.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "nsclc_multi_MS_2.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "nsclc_multi_MS_3.png"
  )

  saveRDS(modscore.demo$harm.object,
          file = test_path("output/modScore/nsclc_multi_so.rds"))

  expect_equal_to_reference(modscore.demo$harm.object, file =
                              test_path("output/modScore/nsclc_multi_so.rds"))

  expected_elements <- c("Likely_CellType",nsclc.multi$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
         failure_message = "modscore results not found")

})

test_that("modScore returns metadata with scores and cell calls for brca", {

  brca <- getModuleScoreParam("brca")

  suppressWarnings(modscore.demo <- do.call(modScore, brca))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "brca_MS_1.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "brca_MS_2.png"
  )

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "brca_MS_3.png"
  )

  saveRDS(modscore.demo$harm.object,
          file = test_path("output/modScore/brca_so.rds"))

  expect_equal_to_reference(modscore.demo$harm.object, file =
                              test_path("output/modScore/brca_so.rds"))

  expected_elements <- c("Likely_CellType",brca$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
         failure_message = "modscore results not found")

})

# Error Testings
test_that("modScore detects when no genes are found in the data", {

  chariou <- getModuleScoreParam("chariou")

  expect_error(suppressWarnings(modscore.demo <- modScore(
    object = chariou$object,
    samples.subset = chariou$samples.subset,
    sample.to.display = chariou$sample.to.display,
    marker.table = apply(chariou$marker.table,2, function(x) toupper(x)),
    celltypes = chariou$celltypes,
    general.class = chariou$general.class,
    lvl.df = chariou$lvl.df,
    nbins = 10), "No genes from list was found in data"))

})

test_that("modScore runs hierarchical classification", {
  
  tec <- getModuleScoreParam("tec")
  
  modscore.demo <- modScore(
    object = tec$object,
    samples.subset = tec$samples.subset,
    sample.to.display = tec$sample.to.display,
    marker.table = tec$marker.table,
    celltypes = c("CD8_T","CD4_T","Tregs"),
    multi.lvl = TRUE,
    general.class = c("CD8_T","CD4_T"),
    lvl.df = tec$lvl.df,
    threshold = c(0.1,0.1,0.1),
    nbins = 10)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[4]]),
    "tec_MS_hierarchical_1.png"
  )
  
  saveRDS(modscore.demo$harm.object,
          file = test_path("output/modScore/TEC_so.rds"))
  
  expect_equal_to_reference(modscore.demo$harm.object,
                            file = test_path("output/modScore/TEC_so.rds"))
  
  expected_elements <- c("Likely_CellType",tec$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("modScore detects when threshold number does not match number of
          cells to analyze", {

  chariou <- getModuleScoreParam("chariou")

  expect_error(suppressWarnings(modscore.demo <- modScore(
    object = chariou$object,
    samples.subset = chariou$samples.subset,
    sample.to.display = chariou$sample.to.display,
    marker.table = chariou$marker.table,
    celltypes = chariou$celltypes,
    general.class = chariou$general.class,
    lvl.df = chariou$lvl.df,
    threshold = rep(0.1,5),
    nbins = 10),
    "Threshold length does not match # celltypes to analyze"))

})
