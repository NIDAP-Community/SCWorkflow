quietPlotMetadata <- function(params) {
  invisible(capture.output(
    output <- tryCatch(do.call(plotMetadata, params), error = identity)
  ))
  if (inherits(output, "error")) {
    stop(conditionMessage(output), call. = FALSE)
  }
  output
}

test_that("plotMetadata prints metadata options and errors when metadata.to.plot is empty", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$samples.to.include <- "c()"
  chariou.data$metadata.to.plot <- "c()"

  captured.output <- capture.output(
    error <- expect_error(
      do.call(plotMetadata, chariou.data),
      "metadata.to.plot must include at least one metadata column"
    )
  )
  log.lines <- c(captured.output, conditionMessage(error))

  expect_true(any(grepl("Found orig.ident in column 1 of object metadata.", captured.output, fixed = TRUE)))
  expect_true(any(grepl("No samples specified. Using all samples...", captured.output, fixed = TRUE)))
  expect_equal(sum(captured.output == "Possible metadata columns to select:"), 1)
  expect_true(any(captured.output == "  - SCT_snn_res.2.4"))
  expect_false(any(captured.output == "  - SCT_snn_res_2_4"))
  expect_equal(
    tail(log.lines, 1),
    "metadata.to.plot must include at least one metadata column. See Possible metadata columns to select above."
  )
})

test_that("plotMetadata errors when object has neither orig.ident nor orig_ident", {
  chariou.data <- getParamPM("Chariou")
  colnames(chariou.data$object@meta.data)[colnames(chariou.data$object@meta.data) == "orig.ident"] <-
    "sample_id"
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'

  expect_error(
    do.call(plotMetadata, chariou.data),
    "Found neither orig.ident nor orig_ident in object metadata. Please provide an object with one of these metadata column names.",
    fixed = TRUE
  )
})

test_that("plotMetadata errors when requested samples are not in the object", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$samples.to.include <- 'c("unknown_sample")'

  invisible(capture.output(
    expect_error(
      do.call(plotMetadata, chariou.data),
      "samples.to.include contains sample names not found in the object: unknown_sample.",
      fixed = TRUE
    )
  ))
})

test_that("plotMetadata errors on metadata columns not found without renaming", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res_2_4")'

  captured.output <- capture.output(
    error <- expect_error(
      do.call(plotMetadata, chariou.data),
      "metadata.to.plot contains metadata columns not found in the object"
    )
  )
  log.lines <- c(captured.output, conditionMessage(error))

  expect_equal(sum(captured.output == "Possible metadata columns to select:"), 1)
  expect_true(any(captured.output == "  - SCT_snn_res.2.4"))
  expect_false(any(captured.output == "  - SCT_snn_res_2_4"))
  expect_equal(
    tail(log.lines, 1),
    "metadata.to.plot contains metadata columns not found in the object: SCT_snn_res_2_4. See Possible metadata columns to select above."
  )
})

test_that("plotMetadata errors when summary columns are not in the object", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$columns.to.summarize <- 'c("missing_summary_column")'

  invisible(capture.output(
    expect_error(
      do.call(plotMetadata, chariou.data),
      "columns.to.summarize contains metadata columns not found in the object: missing_summary_column.",
      fixed = TRUE
    )
  ))
})

test_that("plotMetadata errors on unsupported reduction types", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$reduction.type <- "invalid_reduction"

  expect_error(
    do.call(plotMetadata, chariou.data),
    "reduction.type must be one of: tsne, umap, pca.",
    fixed = TRUE
  )
})

test_that("plotMetadata validates summarization.cut.off", {
  invalid.cutoffs <- list(0, -1, 1.5, NA_real_, Inf, "five")

  for (cutoff in invalid.cutoffs) {
    chariou.data <- getParamPM("Chariou")
    chariou.data$summarization.cut.off <- cutoff

    expect_error(
      do.call(plotMetadata, chariou.data),
      "summarization.cut.off must be a single positive whole number.",
      fixed = TRUE
    )
  }
})

test_that("plotMetadata requires summary cutoffs below unique-value counts", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$object@meta.data$summary_group <- rep(
    c("group_a", "group_b"),
    length.out = nrow(chariou.data$object@meta.data)
  )
  chariou.data$columns.to.summarize <- 'c("summary_group")'
  chariou.data$summarization.cut.off <- 2

  invisible(capture.output(
    expect_error(
      do.call(plotMetadata, chariou.data),
      "summarization.cut.off (2) must be less than the number of unique values (2) in columns.to.summarize column 'summary_group'.",
      fixed = TRUE
    )
  ))
})

test_that("plotMetadata errors when selected metadata column contains NA values", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$object@meta.data$metadata_with_na <- "present"
  chariou.data$object@meta.data$metadata_with_na[1] <- NA
  chariou.data$metadata.to.plot <- 'c("metadata_with_na")'

  captured.output <- capture.output(
    error <- expect_error(
      do.call(plotMetadata, chariou.data),
      "End of error message."
    )
  )

  expect_equal(conditionMessage(error), "End of error message.")
  expect_true(any(grepl("ERROR: Metadata column appears to contain NA values", captured.output)))
  expect_true(any(grepl("Please review your selected metadata column", captured.output, fixed = TRUE)))
  expect_true(any(grepl("Below are valid metadata to select for this plot:", captured.output, fixed = TRUE)))
})

test_that("plotMetadata preserves exact metadata column names", {
  chariou.data <- getParamPM("Chariou")
  original.metadata.columns <- colnames(chariou.data$object@meta.data)
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'

  captured.output <- capture.output(
    output <- do.call(plotMetadata, chariou.data)
  )

  expect_length(output$plots, 1)
  expect_identical(colnames(output$object@meta.data), original.metadata.columns)
  expect_true(any(captured.output == "  - SCT_snn_res.2.4"))
  expect_false(any(captured.output == "  - SCT_snn_res_2_4"))
})

test_that("plotMetadata accepts orig_ident without renaming metadata columns", {
  chariou.data <- getParamPM("Chariou")
  colnames(chariou.data$object@meta.data)[colnames(chariou.data$object@meta.data) == "orig.ident"] <-
    "orig_ident"
  original.metadata.columns <- colnames(chariou.data$object@meta.data)
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'

  captured.output <- capture.output(
    output <- do.call(plotMetadata, chariou.data)
  )

  expect_length(output$plots, 1)
  expect_identical(colnames(output$object@meta.data), original.metadata.columns)
  expect_true("orig_ident" %in% colnames(output$object@meta.data))
  expect_false("orig.ident" %in% colnames(output$object@meta.data))
  expect_true(any(grepl('^\\[1\\] "Found orig_ident in object metadata[.]"$', captured.output)))
})

test_that("plotMetadata adds labels only when show.labels is TRUE", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'

  chariou.data$show.labels <- FALSE
  output.without.labels <- quietPlotMetadata(chariou.data)
  layer.geoms.without.labels <- vapply(
    output.without.labels$plots[[1]]$layers,
    function(layer) class(layer$geom)[1],
    character(1)
  )
  expect_false("GeomLabelRepel" %in% layer.geoms.without.labels)

  chariou.data$show.labels <- TRUE
  output.with.labels <- quietPlotMetadata(chariou.data)
  layer.geoms.with.labels <- vapply(
    output.with.labels$plots[[1]]$layers,
    function(layer) class(layer$geom)[1],
    character(1)
  )
  expect_true("GeomLabelRepel" %in% layer.geoms.with.labels)
})

test_that("plotMetadata applies legend text size and position", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'
  chariou.data$legend.text.size <- 1.5
  chariou.data$legend.position <- "bottom"

  output <- quietPlotMetadata(chariou.data)
  plot <- output$plots[[1]]

  expect_identical(plot$theme$legend.position, "bottom")
  expect_equal(plot$theme$legend.text$size, ggplot2::rel(1.5))
})

test_that("Test Plot Metadata using TEC (Mouse) dataset", {
  tec.data <- getParamPM("TEC")
  output <- quietPlotMetadata(tec.data)

  ggsave("output/TEC_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet1.png")
  ggsave("output/TEC_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet2.png")
  ggsave("output/TEC_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet3.png")
  ggsave("output/TEC_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet4.png")
  ggsave("output/TEC_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet5.png")
  ggsave("output/TEC_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet6.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using Chariou (Mouse) dataset", {
  chariou.data <- getParamPM("Chariou")
  output <- quietPlotMetadata(chariou.data)

  ggsave("output/Chariou_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output", "Chariou_plotmet1.png")
  ggsave("output/Chariou_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output", "Chariou_plotmet2.png")
  ggsave("output/Chariou_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output", "Chariou_plotmet3.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using BRCA (Human) dataset", {
  brca.data <- getParamPM("BRCA")
  output <- quietPlotMetadata(brca.data)

  ggsave("output/BRCA_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output", "BRCA_plotmet1.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output", "BRCA_plotmet2.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output", "BRCA_plotmet3.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output", "BRCA_plotmet4.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output", "BRCA_plotmet5.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output", "BRCA_plotmet6.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using NSCLCmulti (Human) dataset", {
  nsclc.multi.data <- getParamPM("nsclc-multi")
  output <- quietPlotMetadata(nsclc.multi.data)

  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output", "NSCLCmulti_plotmet1.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output", "NSCLCmulti_plotmet2.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output", "NSCLCmulti_plotmet3.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output", "NSCLCmulti_plotmet4.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output", "NSCLCmulti_plotmet5.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output", "NSCLCmulti_plotmet6.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using PBMCsingle (Human) dataset", {
  pbmc.single.data <- getParamPM("pbmc-single")
  output <- quietPlotMetadata(pbmc.single.data)

  ggsave("output/PBMCsingle_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output", "PBMCsingle_plotmet1.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output", "PBMCsingle_plotmet2.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output", "PBMCsingle_plotmet3.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output", "PBMCsingle_plotmet4.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output", "PBMCsingle_plotmet5.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output", "PBMCsingle_plotmet6.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using TEC (Mouse) dataset; UMAP", {
  tec.data <- getParamPM("TEC")
  tec.data$reduction.type <- "umap"
  output <- quietPlotMetadata(tec.data)

  ggsave("output/TEC_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet1.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet2.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet3.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet4.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet5.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output", "TEC_plotmet6.umap.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})
