# Test 1: Standard parameters - TEC dataset
test_that("compareCellPopulations returns correct structure with TEC data", {
  params <- getParamCCP("TEC")
  result <- do.call(compareCellPopulations, params)
  
  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("Plots", "Table"))
  expect_named(result$Plots, c("Barplot", "Boxplot"))
  
  # Check plot types
  expect_s3_class(result$Plots$Barplot, "gg")
  expect_s3_class(result$Plots$Boxplot, "gg")
  
  # Check table structure
  expect_true(is.data.frame(result$Table))
  expect_true("Clusters" %in% colnames(result$Table))
  
  # Snapshot tests for plots and table
  skip_on_ci()
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Barplot),
    "TEC_Standard_Barplot.png"
  )
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Boxplot),
    "TEC_Standard_Boxplot.png"
  )
  expect_snapshot_file(
    .saveCCPTable(result$Table),
    "TEC_Standard_Table.rds"
  )
})

# Test 2: Standard parameters - Chariou dataset
test_that("compareCellPopulations works with Chariou data", {
  params <- getParamCCP("Chariou")
  result <- do.call(compareCellPopulations, params)
  
  expect_type(result, "list")
  expect_s3_class(result$Plots$Barplot, "gg")
  expect_s3_class(result$Plots$Boxplot, "gg")
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Barplot),
    "Chariou_Standard_Barplot.png"
  )
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Boxplot),
    "Chariou_Standard_Boxplot.png"
  )
})

# Test 3: Standard parameters - PBMC dataset with annotated cell types
test_that("compareCellPopulations works with PBMC annotated cell types", {
  params <- getParamCCP("PBMC")
  result <- do.call(compareCellPopulations, params)
  
  expect_type(result, "list")
  expect_s3_class(result$Plots$Barplot, "gg")
  expect_s3_class(result$Plots$Boxplot, "gg")
  
  skip_on_ci()
  # Note: PBMC has only one sample, which creates issues with alluvial flow visualization
  # Skip barplot snapshot for single-sample dataset
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Boxplot),
    "PBMC_Standard_Boxplot.png"
  )
})

# Test 4: Standard parameters - NSCLC dataset
test_that("compareCellPopulations works with NSCLC multi data", {
  params <- getParamCCP("NSCLC")
  result <- do.call(compareCellPopulations, params)
  
  expect_type(result, "list")
  expect_s3_class(result$Plots$Barplot, "gg")
  expect_s3_class(result$Plots$Boxplot, "gg")
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Barplot),
    "NSCLC_Standard_Barplot.png"
  )
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Boxplot),
    "NSCLC_Standard_Boxplot.png"
  )
})

# Test 5: Standard parameters - BRCA dataset
test_that("compareCellPopulations works with BRCA data", {
  params <- getParamCCP("BRCA")
  result <- do.call(compareCellPopulations, params)
  
  expect_type(result, "list")
  expect_s3_class(result$Plots$Barplot, "gg")
  expect_s3_class(result$Plots$Boxplot, "gg")
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Barplot),
    "BRCA_Standard_Barplot.png"
  )
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Boxplot),
    "BRCA_Standard_Boxplot.png"
  )
})

# Test 6: Counts type parameter - TEC dataset
test_that("compareCellPopulations works with Counts type on TEC data", {
  params <- getParamCCP("TEC")
  params$counts.type <- "Counts"
  
  result <- do.call(compareCellPopulations, params)
  
  # Check result structure
  expect_type(result, "list")
  expect_s3_class(result$Plots$Barplot, "gg")
  expect_s3_class(result$Plots$Boxplot, "gg")
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Barplot),
    "TEC_Counts_Barplot.png"
  )
})

# Test 7: Custom group order - Chariou dataset
test_that("compareCellPopulations handles custom group order on Chariou", {
  params <- getParamCCP("Chariou")
  params$group.order <- c("1", "0")  # Status levels
  
  result <- do.call(compareCellPopulations, params)
  
  expect_type(result, "list")
  expect_s3_class(result$Plots$Barplot, "gg")
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Barplot),
    "Chariou_CustomOrder_Barplot.png"
  )
})

# Test 8: Custom wrap columns - PBMC dataset
test_that("compareCellPopulations handles custom wrap.ncols on PBMC", {
  params <- getParamCCP("PBMC")
  params$wrap.ncols <- 3  # Change from default 5 to 3 columns
  
  result <- do.call(compareCellPopulations, params)
  
  expect_type(result, "list")
  expect_s3_class(result$Plots$Barplot, "gg")
  expect_s3_class(result$Plots$Boxplot, "gg")
  
  skip_on_ci()
  expect_snapshot_file(
    .drawCCPFig(result$Plots$Boxplot),
    "PBMC_CustomWrap_Boxplot.png"
  )
})

# Test 9: Input validation - non-Seurat object
test_that("compareCellPopulations validates input object", {
  expect_error(
    compareCellPopulations(
      object = list(),
      metadata.table = data.frame(),
      annotation.column = "cell_type",
      group.column = "treatment"
    ),
    "must be a Seurat object"
  )
})

# Test 10: Missing column validation - TEC dataset
test_that("compareCellPopulations validates metadata columns on TEC", {
  params <- getParamCCP("TEC")
  params$annotation.column <- "nonexistent_column"
  
  expect_error(
    do.call(compareCellPopulations, params),
    "missing from metadata"
  )
})

# Test 11: Invalid counts.type parameter
test_that("compareCellPopulations validates counts.type parameter", {
  params <- getParamCCP("TEC")
  params$counts.type <- "Invalid"
  
  expect_error(
    do.call(compareCellPopulations, params),
    "must be either 'Frequency' or 'Counts'"
  )
})

# Test 12: Table output validation - BRCA dataset
test_that("compareCellPopulations table contains expected columns on BRCA", {
  params <- getParamCCP("BRCA")
  result <- do.call(compareCellPopulations, params)
  
  # Check for _CellCounts and _Percent suffixed columns
  expect_true(any(grepl("_CellCounts$", colnames(result$Table))))
  expect_true(any(grepl("_Percent$", colnames(result$Table))))
})

