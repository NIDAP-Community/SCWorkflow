test_that("Harmony returns seurat object with adjusted embeddings for
          TEC data", {
  skip_on_ci()

  tec = getHarmonyParam("TEC")

  object.harmonized = do.call(harmonyBatchCorrect, tec)
   
  # expect_snapshot_file(
  #   .drawHarmonyFig(object.harmonized$adj.tsne),
  #   "tec_harm.png"
  # )
  
  expect_equal(mean(object.harmonized[["Harmony"]]["Lum",]),
               0.06279209, tolerance = 5e-1)

  # expected_elements = c("adj.object","adj.tsne")
  # expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for Chariou
          data", {
  skip_on_ci()
  chariou = getHarmonyParam("Chariou")

  object.harmonized = do.call(harmonyBatchCorrect, chariou)


  # expect_snapshot_file(
  #   .drawHarmonyFig(object.harmonized$adj.tsne),
  #   "char_harm.png"
  # )

  expect_equal(mean(object.harmonized[["Harmony"]]["Ccl8",]),
               0.2046126, tolerance = 5e-1)

  #expected_elements = c("adj.object","adj.tsne")
  #expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          pbmc_single data", {
  skip_on_ci()
  pbmc.single = getHarmonyParam("pbmc_single")

  object.harmonized = do.call(harmonyBatchCorrect, pbmc.single)

  
  # expect_snapshot_file(
  #   .drawHarmonyFig(object.harmonized$adj.tsne),
  #   "pbmc_single_harm.png"
  # )

  expect_equal(mean(object.harmonized[["Harmony"]]["CLU",]),
               -0.002671558, tolerance = 5e-1)

  # expected_elements = c("adj.object","adj.tsne")
  # expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          nsclc_multi data", {
  skip_on_ci()
  nsclc.multi = getHarmonyParam("nsclc_multi")

  object.harmonized = do.call(harmonyBatchCorrect, nsclc.multi)

  # expect_snapshot_file(
  #   .drawHarmonyFig(object.harmonized$adj.tsne),
  #   "nsclc_multi_harm.png"
  # )

  expect_equal(mean(object.harmonized[["Harmony"]]["LCN2",]),
               0.1265227, tolerance = 5e-1)

  # expected_elements = c("adj.object","adj.tsne")
  # expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          BRCA data", {
  skip_on_ci()
          
  brca = getHarmonyParam("BRCA")

  object.harmonized = do.call(harmonyBatchCorrect, brca)

  # expect_snapshot_file(
  #   .drawHarmonyFig(object.harmonized$adj.tsne),
  #   "BRCA_harm.png"
  # )

  expect_equal(mean(object.harmonized[["Harmony"]]["SCGB2A2",]),
               -0.1078249, tolerance = 5e-1)

  # expected_elements = c("adj.object","adj.tsne")
  # expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony provides warning when genes are not found in the data", {
  skip_on_ci()
  tec = getHarmonyParam("TEC")

  expect_warning(harmonyBatchCorrect(
    object = tec$object,
    nvar = tec$nvar,
    genes.to.add = c("wrong_gene","wrong_gene2"),
    group.by.var = tec$group.by.var),
    "specified genes were not found and therefore cannot be added")

})

test_that("Harmony stops when variable features to subset by exceeds number of
          genes in the data", {
  skip_on_ci()
  tec = getHarmonyParam("TEC")

  expect_error(harmonyBatchCorrect(
    object = tec$object,
    nvar = 100000,
    genes.to.add = NULL,
    group.by.var = tec$group.by.var),
    "nvar exceed total number of variable genes in the data")

})