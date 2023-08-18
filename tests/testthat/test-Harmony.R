test_that("Harmony returns seurat object with adjusted embeddings for
          TEC data", {

  tec = getHarmonyParam("TEC")

  object.harmonized = do.call(harmonyBatchCorrect, tec)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$harm.figures),
    "tec_harm.png"
  )

  saveRDS(object.harmonized$harm.object,
          file = test_path("output/harmony/TEC_so.rds"))

  expect_equal_to_reference(object.harmonized$harm.object,
                            file = test_path("output/harmony/TEC_so.rds"))

  expected_elements = c("harm.object","harm.figures")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for Chariou
          data", {

  chariou = getHarmonyParam("Chariou")

  object.harmonized = do.call(harmonyBatchCorrect, chariou)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$harm.figures),
    "chariou_harm.png"
  )

  saveRDS(object.harmonized$harm.object,
          file = test_path("output/harmony/chariou_so.rds"))

  expect_equal_to_reference(object.harmonized$harm.object,
                            file = test_path("output/harmony/chariou_so.rds"))

  expected_elements = c("harm.object","harm.figures")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          pbmc_single data", {

  pbmc.single = getHarmonyParam("pbmc_single")

  object.harmonized = do.call(harmonyBatchCorrect, pbmc.single)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$harm.figures),
    "pbmc_single_harm.png"
  )

  saveRDS(object.harmonized$harm.object,
          file = test_path("output/harmony/pbmc_single_so.rds"))

  expect_equal_to_reference(object.harmonized$harm.object,
                            file = test_path("output/harmony/pbmc_single_so.rds"))

  expected_elements = c("harm.object","harm.figures")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          nsclc_multi data", {

  nsclc_multi = getHarmonyParam("nsclc_multi")

  object.harmonized = do.call(harmonyBatchCorrect, nsclc_multi)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$harm.figures),
    "nsclc_multi_harm.png"
  )

  saveRDS(object.harmonized$harm.object,
          file = test_path("output/harmony/nsclc_multi_so.rds"))

  expect_equal_to_reference(object.harmonized$harm.object, file = 
                              test_path("output/harmony/nsclc_multi_so.rds"))

  expected_elements = c("harm.object","harm.figures")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          BRCA data", {

  brca = getHarmonyParam("BRCA")

  object.harmonized = do.call(harmonyBatchCorrect, brca)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$harm.figures),
    "brca_harm.png"
  )

  saveRDS(object.harmonized$harm.object,
          file = test_path("output/harmony/brca_so.rds"))

  expect_equal_to_reference(object.harmonized$harm.object,
                            file = test_path("output/harmony/brca_so.rds"))

  expected_elements = c("harm.object","harm.figures")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony provides warning when genes are not found in the data", {

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

  tec = getHarmonyParam("TEC")

  expect_error(harmonyBatchCorrect(
    object = tec$object,
    nvar = 100000,
    genes.to.add = NULL,
    group.by.var = tec$group.by.var),
    "nvar exceed total number of variable genes in the data")

})