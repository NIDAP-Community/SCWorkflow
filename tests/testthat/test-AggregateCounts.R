test_that("aggregateCounts returns pseudobulk with sanitized group columns", {
  obj <- selectCRObject("BRCA")

  res <- aggregateCounts(
    object = obj,
    var.group = "orig.ident",
    slot = "data"
  )

  expect_true(is.data.frame(res))
  expect_true("Gene" %in% colnames(res))
  # Row names are not relied upon; ensure Gene column exists instead

  meta_groups <- unique(as.character(obj$orig.ident))
  expected_cols <- gsub("\\W", "_", meta_groups)

  expect_setequal(setdiff(colnames(res), "Gene"), expected_cols)
})

test_that("aggregateCounts warns for singleton groups", {
  obj <- selectCRObject("BRCA")

  # Force a singleton group in orig.ident
  orig <- as.character(obj$orig.ident)
  orig[1] <- paste0(orig[1], "_ONE")
  obj$orig.ident <- factor(orig)

  expect_warning(
    res <- aggregateCounts(
      object = obj,
      var.group = "orig.ident",
      slot = "data"
    ),
    "Some groups have only 1 cell"
  )

  # Ensure the new singleton group column exists (sanitized)
  expected_col <- gsub("\\W", "_", orig[1])
  expect_true(expected_col %in% colnames(res))
})

test_that("aggregateCounts errors for non-categorical var.group", {
  obj <- selectCRObject("BRCA")
  obj$numeric_group <- seq_len(ncol(obj))

  expect_error(
    aggregateCounts(
      object = obj,
      var.group = "numeric_group",
      slot = "data"
    ),
    "All columns in var.group must be factors or characters"
  )
})
