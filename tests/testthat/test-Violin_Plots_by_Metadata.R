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

test_that("Violin plot works for NSCLC_Single data", {
  NSCLC_Single = selectViolin("NSCLC_Single")

  violin_test = do.call(violinPlot, NSCLC_Single)

  skip_on_ci()
  expect_snapshot_file(
    .drawViolin(violin_test$plots),
    "NSCLC_Single_violin.png"
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
#   NSCLC_Single <- selectViolin("NSCLC_Single")
# 
#   expect_error(
#     violinPlot(
#       object = NSCLC_Single$object,
#       group.by = NSCLC_Single$group.by,
#       group.subset = NSCLC_Single$group.subset,
#       genes.of.interest = paste("jibberish", 1:5, sep =
#                                   "_")
#     ),
#     "No query genes were found in the dataset."
#   )
# 
# })

# test_that("Violin plot stops when ident of interest is not found in seurat", {
#   NSCLC_Single <- selectViolin("NSCLC_Single")
# 
#   expect_error(
#     violinPlot(
#       object = NSCLC_Single$object,
#       group.by = "jibberish",
#       group.subset = NSCLC_Single$group.subset,
#       genes.of.interest = NSCLC_Single$genes.of.interest
#     ),
#     "Unable to find ident of interest in metadata."
#   )
# 
# })

# test_that("Violin plot stops when user attempts to rename group.by as
#           Gene, Expression, or Scaled",
#           {
#             NSCLC_Single <- selectViolin("NSCLC_Single")
# 
#             expect_error(
#               violinPlot(
#                 object = NSCLC_Single$object,
#                 group.by = NSCLC_Single$group.by,
#                 group.subset = NSCLC_Single$group.subset,
#                 genes.of.interest = NSCLC_Single$genes.of.interest,
#                 rename.ident = "Gene"
#               ),
#               "New ident name cannot be one of Gene, Expression, or scaled."
#             )
# 
#           })
