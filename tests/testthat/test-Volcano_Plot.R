# Tests for Volcano_Plot function

library(testthat)

# Helper function to create test data
create_test_deg_data <- function(n_genes = 100) {
  set.seed(123)
  data.frame(
    Gene = paste0("Gene", 1:n_genes),
    `WT-KO_logFC` = rnorm(n_genes, 0, 2),
    `WT-KO_pval` = runif(n_genes, 0, 0.1),
    `Treated-Control_logFC` = rnorm(n_genes, 0, 1.5),
    `Treated-Control_pval` = runif(n_genes, 0, 0.05),
    check.names = FALSE
  )
}

test_that("Volcano_Plot runs without errors on valid input", {
  deg_data <- create_test_deg_data()
  
  expect_no_error(
    result <- Volcano_Plot(
      DEGAnalysis = deg_data,
      label.col = "Gene",
      sig.col = "WT-KO_pval",
      lfc.col = "WT-KO_logFC",
      pCutoff = 0.05,
      FCcutoff = 1.0
    )
  )
})

test_that("Volcano_Plot returns correct structure", {
  deg_data <- create_test_deg_data()
  
  result <- Volcano_Plot(
    DEGAnalysis = deg_data,
    label.col = "Gene",
    sig.col = "WT-KO_pval",
    lfc.col = "WT-KO_logFC",
    pCutoff = 0.05,
    FCcutoff = 1.0
  )
  
  # Check that result is a list
  expect_type(result, "list")
  
  # Check that 'data' element exists
  expect_true("data" %in% names(result))
  
  # Check that data is a data frame
  expect_s3_class(result$data, "data.frame")
  
  # Check that rank column was added
  expect_true(any(grepl("_rank$", colnames(result$data))))
  
  # Check that static plot exists
  expect_true(any(grepl("_static$", names(result))))
  
  # Check that interactive plot exists
  expect_true(any(grepl("_interactive$", names(result))))
  
  # Check that static plot is a ggplot object
  static_plot_name <- names(result)[grepl("_static$", names(result))][1]
  expect_s3_class(result[[static_plot_name]], "ggplot")
  
  # Check that interactive plot is a plotly object
  interactive_plot_name <- names(result)[grepl("_interactive$", names(result))][1]
  expect_s3_class(result[[interactive_plot_name]], "plotly")
})

test_that("Volcano_Plot handles multiple comparisons correctly", {
  deg_data <- create_test_deg_data()
  
  result <- Volcano_Plot(
    DEGAnalysis = deg_data,
    label.col = "Gene",
    sig.col = c("WT-KO_pval", "Treated-Control_pval"),
    lfc.col = c("WT-KO_logFC", "Treated-Control_logFC"),
    pCutoff = 0.05,
    FCcutoff = 1.0
  )
  
  # Should have data + 4 plots (2 comparisons × 2 plot types)
  expect_equal(length(result), 5)
  
  # Check that both comparison names appear
  expect_true(any(grepl("WT-KO", names(result))))
  expect_true(any(grepl("Treated-Control", names(result))))
  
  # Check that we have both static and interactive for each comparison
  expect_equal(sum(grepl("_static$", names(result))), 2)
  expect_equal(sum(grepl("_interactive$", names(result))), 2)
})

test_that("Volcano_Plot handles edge cases", {
  # Test with p-values = 0
  deg_data <- create_test_deg_data(n_genes = 50)
  deg_data$`WT-KO_pval`[1:5] <- 0  # Set some p-values to exactly 0
  
  expect_no_error(
    result <- Volcano_Plot(
      DEGAnalysis = deg_data,
      label.col = "Gene",
      sig.col = "WT-KO_pval",
      lfc.col = "WT-KO_logFC",
      pCutoff = 0.05,
      FCcutoff = 1.0
    )
  )
  
  # Check that the function completed and returned data
  expect_true("data" %in% names(result))
  expect_equal(nrow(result$data), 50)
  
  # Test with NA values
  deg_data_na <- create_test_deg_data(n_genes = 50)
  deg_data_na$`WT-KO_logFC`[1:3] <- NA
  deg_data_na$`WT-KO_pval`[4:6] <- NA
  
  expect_no_error(
    result_na <- Volcano_Plot(
      DEGAnalysis = deg_data_na,
      label.col = "Gene",
      sig.col = "WT-KO_pval",
      lfc.col = "WT-KO_logFC",
      pCutoff = 0.05,
      FCcutoff = 1.0
    )
  )
  
  # Check that the function completed and returned valid data
  # (NA values are handled internally for plotting but preserved in returned data)
  expect_true("data" %in% names(result_na))
  expect_equal(nrow(result_na$data), 50)
  expect_s3_class(result_na$`WT-KO_static`, "ggplot")
  expect_s3_class(result_na$`WT-KO_interactive`, "plotly")
})