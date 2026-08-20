test_that("Violin plot works for TEC data", {
  tec.data = selectViolin("TEC")
  
  violin_test = do.call(violinPlot, tec.data)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test$plots),
    "tec_violin.png"
  )
  
  expect_setequal(any(grepl("ggplot", class(violin_test$plots))), TRUE)
  
})

test_that("Violin plot works for Chariou data", {
  chariou.data = selectViolin("Chariou")

  violin_test = do.call(violinPlot, chariou.data)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test$plots),
    "chariou_violin.png"
  )

  
  expect_setequal(any(grepl("ggplot", class(violin_test$plots))), TRUE)

})

# test_that("Violin plot works for Chariou.allgroup data", {
#   chariou.allgroup.data = selectViolin("Chariou.allgroups")
# 
#   violin_test = do.call(violinPlot, chariou.allgroup.data)
# 
#   skip_on_ci()
#   expect_snapshot_file(
#     .drawViolin(violin_test$plots),
#     "chariou_allgroup_violin.png"
#   )
# 
#   
#   expect_setequal(any(grepl("ggplot", class(violin_test$plots))), TRUE)
# 
# })
# 
# test_that("Violin plot works for Chariou.subgroup data", {
#   chariou.subgroup.data = selectViolin("Chariou.subgroup")
# 
#   violin_test = do.call(violinPlot, chariou.subgroup.data)
# 
#   skip_on_ci()
#   expect_snapshot_file(
#     .drawViolin(violin_test$plots),
#     "chariou_subgroup_violin.png"
#   )
# 
#   
#   expect_setequal(any(grepl("ggplot", class(violin_test$plots))), TRUE)
# 
# })

test_that("Violin plot works for pbmc.single data", {
  pbmc.single = selectViolin("pbmc.single")

  violin_test = do.call(violinPlot, pbmc.single)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test$plots),
    "pbmc_single_violin.png"
  )

  
  expect_setequal(any(grepl("ggplot", class(violin_test$plots))), TRUE)

})

test_that("Violin plot works for nsclc.multi data", {
  nsclc.multi = selectViolin("nsclc.multi")

  violin_test = do.call(violinPlot, nsclc.multi)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test$plots),
    "nsclc_multi_violin.png"
  )

  
  expect_setequal(any(grepl("ggplot", class(violin_test$plots))), TRUE)

})

test_that("Violin plot works for brca data", {
  brca = selectViolin("brca")

  violin_test = do.call(violinPlot, brca)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test$plots),
    "brca_violin.png"
  )

  
  expect_setequal(any(grepl("ggplot", class(violin_test$plots))), TRUE)

})

## Check code detects warnings and errors ##

# test_that("Violin plot stops when no query genes are found in the data", {
#   pbmc.single <- selectViolin("pbmc.single")
# 
#   expect_error(
#     violinPlot(
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
#     violinPlot(
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
#               violinPlot(
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
