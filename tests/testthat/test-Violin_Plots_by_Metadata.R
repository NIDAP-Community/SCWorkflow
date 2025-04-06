test_that("Violin plot works for TEC data", {
  tec.data = selectViolin("TEC")
  
  violin_test = do.call(violinPlot_mod, tec.data)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test),
    "tec_violin.png"
  )
  
  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for Chariou data", {
  chariou.data = selectViolin("Chariou")

  violin_test = do.call(violinPlot_mod, chariou.data)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test),
    "chariou_violin.png"
  )

  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)

})

# test_that("Violin plot works for Chariou.allgroup data", {
#   chariou.allgroup.data = selectViolin("Chariou.allgroups")
# 
#   violin_test = do.call(violinPlot_mod, chariou.allgroup.data)
# 
#   skip_on_ci()
#   expect_snapshot_file(
#     .drawViolin(violin_test),
#     "chariou_allgroup_violin.png"
#   )
# 
#   expected_elements = c("gg", "ggplot")
#   expect_setequal(class(violin_test), expected_elements)
# 
# })
# 
# test_that("Violin plot works for Chariou.subgroup data", {
#   chariou.subgroup.data = selectViolin("Chariou.subgroup")
# 
#   violin_test = do.call(violinPlot_mod, chariou.subgroup.data)
# 
#   skip_on_ci()
#   expect_snapshot_file(
#     .drawViolin(violin_test),
#     "chariou_subgroup_violin.png"
#   )
# 
#   expected_elements = c("gg", "ggplot")
#   expect_setequal(class(violin_test), expected_elements)
# 
# })

test_that("Violin plot works for pbmc.single data", {
  pbmc.single = selectViolin("pbmc.single")

  violin_test = do.call(violinPlot_mod, pbmc.single)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test),
    "pbmc_single_violin.png"
  )

  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)

})

test_that("Violin plot works for nsclc.multi data", {
  nsclc.multi = selectViolin("nsclc.multi")

  violin_test = do.call(violinPlot_mod, nsclc.multi)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test),
    "nsclc_multi_violin.png"
  )

  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)

})

test_that("Violin plot works for brca data", {
  brca = selectViolin("brca")

  violin_test = do.call(violinPlot_mod, brca)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test),
    "brca_violin.png"
  )

  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)

})

## Check code detects warnings and errors ##

# test_that("Violin plot stops when no query genes are found in the data", {
#   pbmc.single <- selectViolin("pbmc.single")
# 
#   expect_error(
#     violinPlot_mod(
#       object = pbmc.single$object,
#       group.by = pbmc.single$group.by,
#       group.subset = pbmc.single$group.subset,
#       genes.of.interest = paste("jibberish", 1:5, sep =
#                                   "_")
#     ),
#     "No query genes were found in the dataset."
#   )
# 
# })

# test_that("Violin plot stops when ident of interest is not found in seurat", {
#   pbmc.single <- selectViolin("pbmc.single")
# 
#   expect_error(
#     violinPlot_mod(
#       object = pbmc.single$object,
#       group.by = "jibberish",
#       group.subset = pbmc.single$group.subset,
#       genes.of.interest = pbmc.single$genes.of.interest
#     ),
#     "Unable to find ident of interest in metadata."
#   )
# 
# })

# test_that("Violin plot stops when user attempts to rename group.by as
#           Gene, Expression, or Scaled",
#           {
#             pbmc.single <- selectViolin("pbmc.single")
# 
#             expect_error(
#               violinPlot_mod(
#                 object = pbmc.single$object,
#                 group.by = pbmc.single$group.by,
#                 group.subset = pbmc.single$group.subset,
#                 genes.of.interest = pbmc.single$genes.of.interest,
#                 rename.ident = "Gene"
#               ),
#               "New ident name cannot be one of Gene, Expression, or scaled."
#             )
# 
#           })
