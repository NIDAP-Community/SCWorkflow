test_that("plotMetadata prints metadata options and errors when metadata.to.plot is empty", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$samples.to.include <- "c()"
  chariou.data$metadata.to.plot <- "c()"

  captured.messages <- character()
  captured.output <- capture.output(
    error <- withCallingHandlers(
      expect_error(
        do.call(plotMetadata, chariou.data),
        "metadata.to.plot must include at least one metadata column"
      ),
      message = function(condition) {
        captured.messages <<- c(captured.messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      }
    )
  )
  log.lines <- c(captured.output, captured.messages, conditionMessage(error))

  expect_true(any(grepl("Found orig.ident in object metadata.", captured.output, fixed = TRUE)))
  expect_true(any(grepl("No samples specified. Using all samples...", captured.output, fixed = TRUE)))
  metadata.options.messages <- captured.messages[
    startsWith(captured.messages, "Possible metadata columns to select:")
  ]
  expect_length(metadata.options.messages, 1)
  expect_match(metadata.options.messages, "Possible metadata columns to select:", fixed = TRUE)
  expect_true(any(grepl("  - SCT_snn_res.2.4", captured.messages, fixed = TRUE)))
  expect_false(any(grepl("  - SCT_snn_res_2_4", captured.messages, fixed = TRUE)))
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

  expect_error(
    do.call(plotMetadata, chariou.data),
    "samples.to.include contains sample names not found in the object: unknown_sample.",
    fixed = TRUE
  )
})

test_that("plotMetadata lists rejected and available samples for invalid selections", {
  for (sample.metadata.column in c("orig.ident", "orig_ident")) {
    input <- getParamPM("Chariou")
    colnames(input$object@meta.data)[colnames(input$object@meta.data) == "orig.ident"] <-
      sample.metadata.column
    available.samples <- unique(input$object@meta.data[[sample.metadata.column]])
    invalid.samples <- c("unknown_sample", "another_missing_sample", "unknown_sample")

    for (selected.samples in list(
      invalid.samples,
      c(as.character(available.samples[[1]]), invalid.samples)
    )) {
      input$samples.to.include <- selected.samples

      error <- expect_error(
        do.call(plotMetadata, input),
        "samples.to.include contains sample names not found in the object: unknown_sample, another_missing_sample.",
        fixed = TRUE
      )
      expect_identical(
        conditionMessage(error),
        paste0(
          "samples.to.include contains sample names not found in the object: unknown_sample, another_missing_sample.\n",
          "Available samples (", sample.metadata.column, "):\n",
          paste0("  - ", available.samples, collapse = "\n")
        )
      )
    }
  }
})

test_that("plotMetadata logs selected samples and preserves the full returned dataset", {
  for (sample.metadata.column in c("orig.ident", "orig_ident")) {
    input <- getParamPM("Chariou")
    colnames(input$object@meta.data)[colnames(input$object@meta.data) == "orig.ident"] <-
      sample.metadata.column
    original.metadata <- input$object@meta.data
    selected.samples <- "CD8dep"
    selected.cells <- rownames(original.metadata)[
      original.metadata[[sample.metadata.column]] %in% selected.samples
    ]
    input$samples.to.include <- selected.samples
    input$metadata.to.plot <- sample.metadata.column

    expect_gt(length(selected.cells), 0)
    expect_lt(length(selected.cells), ncol(input$object))

    captured.messages <- character()
    invisible(capture.output(
      output <- withCallingHandlers(
        do.call(plotMetadata, input),
        message = function(condition) {
          captured.messages <<- c(captured.messages, conditionMessage(condition))
          invokeRestart("muffleMessage")
        }
      )
    ))

    expect_length(output$plots, 1)
    expect_equal(nrow(output$plots[[1]]$data), length(selected.cells))
    expect_setequal(output$plots[[1]]$data$clusid, selected.samples)
    expect_equal(
      output$plots[[1]]$data$umap1,
      unname(Seurat::Embeddings(input$object, "umap")[selected.cells, 1])
    )
    expect_identical(colnames(output$object), colnames(input$object))
    expect_identical(output$object@meta.data, original.metadata)
    expect_identical(
      captured.messages[[1]],
      paste0("Selected samples (", sample.metadata.column, "): CD8dep\n")
    )
    expect_identical(
      captured.messages[[2]],
      sprintf("Plotting %s of %s cells.\n", length(selected.cells), ncol(input$object))
    )
  }
})

test_that("plotMetadata logs all cells when no samples are selected", {
  for (sample.metadata.column in c("orig.ident", "orig_ident")) {
    input <- getParamPM("Chariou")
    colnames(input$object@meta.data)[colnames(input$object@meta.data) == "orig.ident"] <-
      sample.metadata.column
    input$metadata.to.plot <- sample.metadata.column
    total.cells <- ncol(input$object)
    available.samples <- unique(input$object@meta.data[[sample.metadata.column]])

    for (empty.selection in list(NULL, character(), "c()")) {
      input["samples.to.include"] <- list(empty.selection)
      captured.messages <- character()
      captured.output <- capture.output(
        output <- withCallingHandlers(
          do.call(plotMetadata, input),
          message = function(condition) {
            captured.messages <<- c(captured.messages, conditionMessage(condition))
            invokeRestart("muffleMessage")
          }
        )
      )

      expect_true(any(grepl(
        "No samples specified. Using all samples...", captured.output, fixed = TRUE
      )))
      expect_equal(nrow(output$plots[[1]]$data), total.cells)
      expect_setequal(output$plots[[1]]$data$clusid, as.character(available.samples))
      expect_identical(output$object@meta.data, input$object@meta.data)
      expect_identical(
        captured.messages[[1]],
        paste0(
          "Selected samples (", sample.metadata.column, "): ",
          paste(available.samples, collapse = ", "), "\n"
        )
      )
      expect_identical(
        captured.messages[[2]],
        sprintf("Plotting %s of %s cells.\n", total.cells, total.cells)
      )
    }
  }
})

test_that("plotMetadata errors on metadata columns not found without renaming", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res_2_4")'

  captured.messages <- character()
  captured.output <- capture.output(
    error <- withCallingHandlers(
      expect_error(
        do.call(plotMetadata, chariou.data),
        "metadata.to.plot contains metadata columns not found in the object"
      ),
      message = function(condition) {
        captured.messages <<- c(captured.messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      }
    )
  )
  log.lines <- c(captured.output, captured.messages, conditionMessage(error))

  expect_equal(sum(startsWith(captured.messages, "Possible metadata columns to select:")), 1)
  expect_true(any(grepl("  - SCT_snn_res.2.4", captured.messages, fixed = TRUE)))
  expect_false(any(grepl("  - SCT_snn_res_2_4", captured.messages, fixed = TRUE)))
  expect_equal(
    tail(log.lines, 1),
    "metadata.to.plot contains metadata columns not found in the object: SCT_snn_res_2_4. See Possible metadata columns to select above."
  )
})

test_that("plotMetadata errors when summary columns are not in the object", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$columns.to.summarize <- 'c("missing_summary_column")'

  expect_error(
    do.call(plotMetadata, chariou.data),
    "columns.to.summarize contains metadata columns not found in the object: missing_summary_column.",
    fixed = TRUE
  )
})

test_that("plotMetadata errors on unsupported reduction types", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$reduction.type <- "invalid_reduction"

  expect_error(
    do.call(plotMetadata, chariou.data),
    "reduction.type must be one of: tsne, umap, pca.",
    fixed = TRUE
  )

  chariou.data$reduction.type <- c("umap", "tsne")
  expect_error(
    do.call(plotMetadata, chariou.data),
    "reduction.type must be a single character value.",
    fixed = TRUE
  )
})

test_that("plotMetadata reports absent RNA and protein reductions", {
  for (use.cite.seq in c(FALSE, TRUE)) {
    for (reduction.type in c("umap", "tsne", "pca")) {
      input <- getParamPM("Chariou")
      input$reduction.type <- reduction.type
      input$use.cite.seq <- use.cite.seq
      reduction.name <- if (use.cite.seq) paste0("protein_", reduction.type) else reduction.type
      input$object@reductions[[reduction.name]] <- NULL
      available.reductions <- names(input$object@reductions)

      error <- expect_error(
        do.call(plotMetadata, input),
        paste0("Requested reduction '", reduction.name, "' was not found in the Seurat object."),
        fixed = TRUE
      )
      expect_identical(
        conditionMessage(error),
        paste0(
          "Requested reduction '", reduction.name, "' was not found in the Seurat object.\n",
          "Available reductions:\n", paste0("  - ", available.reductions, collapse = "\n")
        )
      )
    }
  }

  input <- getParamPM("Chariou")
  input$object@reductions <- list()
  input$reduction.type <- "umap"
  expect_error(
    do.call(plotMetadata, input),
    "Requested reduction 'umap' was not found in the Seurat object.\nAvailable reductions:\n  - <none>",
    fixed = TRUE
  )
})

test_that("plotMetadata logs resolved reductions after sample counts", {
  for (use.cite.seq in c(FALSE, TRUE)) {
    for (reduction.type in c("umap", "tsne", "pca")) {
      input <- getParamPM("Chariou")
      input$reduction.type <- reduction.type
      input$use.cite.seq <- use.cite.seq
      input$samples.to.include <- "CD8dep"
      input$metadata.to.plot <- "orig.ident"
      reduction.name <- if (use.cite.seq) paste0("protein_", reduction.type) else reduction.type
      if (use.cite.seq) {
        embeddings <- Seurat::Embeddings(input$object, reduction.type)[, 1:2, drop = FALSE]
        colnames(embeddings) <- c("protein_1", "protein_2")
        input$object[[reduction.name]] <- SeuratObject::CreateDimReducObject(
          embeddings = embeddings, key = "protein_",
          assay = SeuratObject::DefaultAssay(input$object)
        )
      }
      selected.cells <- rownames(input$object@meta.data)[input$object@meta.data$orig.ident == "CD8dep"]

      captured.messages <- character()
      invisible(capture.output(
        output <- withCallingHandlers(
          do.call(plotMetadata, input),
          message = function(condition) {
            captured.messages <<- c(captured.messages, conditionMessage(condition))
            invokeRestart("muffleMessage")
          }
        )
      ))
      sample.position <- which(startsWith(captured.messages, "Selected samples (orig.ident):"))
      count.position <- which(startsWith(captured.messages, "Plotting "))
      reduction.position <- which(startsWith(captured.messages, "Running plotMetadata with reduction:"))

      expect_length(reduction.position, 1)
      expect_lt(sample.position, count.position)
      expect_lt(count.position, reduction.position)
      expect_identical(captured.messages[reduction.position], paste0("Running plotMetadata with reduction: ", reduction.name, "\n"))
      expect_equal(nrow(output$plots[[1]]$data), length(selected.cells))
      expect_equal(output$plots[[1]]$data$umap1, unname(Seurat::Embeddings(input$object, reduction.name)[selected.cells, 1]))
      expect_equal(output$plots[[1]]$data$umap2, unname(Seurat::Embeddings(input$object, reduction.name)[selected.cells, 2]))
    }
  }
})

test_that("plotMetadata logs the reduction before attempting to draw plots", {
  input <- getParamPM("Chariou")
  input$reduction.type <- "umap"
  input$metadata.to.plot <- "orig.ident"
  input$object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = Seurat::Embeddings(input$object, "umap")[, 1, drop = FALSE],
    key = "UMAP_", assay = SeuratObject::DefaultAssay(input$object)
  )
  captured.messages <- character()
  invisible(capture.output(
    expect_error(withCallingHandlers(
      do.call(plotMetadata, input),
      message = function(condition) {
        captured.messages <<- c(captured.messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      }
    ))
  ))
  expect_true(any(startsWith(captured.messages, "Plotting ")))
  expect_identical(tail(captured.messages, 1), "Running plotMetadata with reduction: umap\n")
})

test_that("plotMetadata defaults an empty reduction type to umap", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$reduction.type <- character()

  output <- do.call(plotMetadata, chariou.data)

  expect_length(output$plots, 3)
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

  expect_error(
    do.call(plotMetadata, chariou.data),
    "Too many categories were requested for columns.to.summarize column 'summary_group': 2 requested, but 1 categories are possible.",
    fixed = TRUE
  )
})

test_that("plotMetadata validates summary cutoffs against selected samples", {
  chariou.data <- getParamPM("Chariou")
  selected.sample <- "CD8dep"
  selected.cells <- chariou.data$object@meta.data$orig.ident == selected.sample
  chariou.data$object@meta.data$summary_group <- ifelse(
    selected.cells,
    rep(c("group_a", "group_b"), length.out = sum(selected.cells)),
    rep(c("group_c", "group_d"), length.out = sum(!selected.cells))
  )
  chariou.data$samples.to.include <- sprintf('c("%s")', selected.sample)
  chariou.data$columns.to.summarize <- 'c("summary_group")'
  chariou.data$summarization.cut.off <- 2

  expect_error(
    do.call(plotMetadata, chariou.data),
    "Too many categories were requested for columns.to.summarize column 'summary_group': 2 requested, but 1 categories are possible.",
    fixed = TRUE
  )
})

test_that("plotMetadata permits summary cutoffs above 20 when values are available", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$object@meta.data$summary_group <- rep(
    sprintf("group_%02d", 1:22),
    length.out = nrow(chariou.data$object@meta.data)
  )
  chariou.data$columns.to.summarize <- 'c("summary_group")'
  chariou.data$summarization.cut.off <- 21

  output <- do.call(plotMetadata, chariou.data)

  expect_length(output$plots, 4)
})

test_that("plotMetadata logs applied metadata summarization", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$columns.to.summarize <- 'c("SCT_snn_res.2.4")'
  chariou.data$summarization.cut.off <- 5

  captured.output <- capture.output(
    output <- do.call(plotMetadata, chariou.data)
  )

  expect_length(output$plots, 3)
  expect_true(any(captured.output == paste0(
    "[1] \"Summarizing metadata column 'SCT_snn_res.2.4': retaining the 5 most frequent values, labeling non-retained values as Other, and labeling missing values as NA.\""
  )))
})

test_that("plotMetadata retains missing summarized values as grey NA cells", {
  chariou.data <- getParamPM("Chariou")
  chariou.data$object@meta.data$summary_with_missing <- rep(
    c("group_a", "group_b", "group_c", NA_character_),
    length.out = nrow(chariou.data$object@meta.data)
  )
  chariou.data$metadata.to.plot <- 'c("Status")'
  chariou.data$columns.to.summarize <- 'c("summary_with_missing")'
  chariou.data$summarization.cut.off <- 2

  output <- do.call(plotMetadata, chariou.data)

  plot <- output$plots[[2]]
  built.plot <- ggplot2::ggplot_build(plot)
  missing.cells <- is.na(plot$data$clusid)
  expect_true(any(missing.cells))
  expect_equal(nrow(plot$data), ncol(chariou.data$object))
  expect_true(all(built.plot$data[[1]]$colour[missing.cells] == "grey50"))
  expect_true(anyNA(built.plot$plot$scales$get_scales("colour")$get_breaks()))
  expect_identical(output$object@meta.data, chariou.data$object@meta.data)
})

test_that("plotMetadata retains categorical NA cells in grey with an NA legend entry", {
  input <- getParamPM("Chariou")
  input$metadata.to.plot <- "metadata_with_na"
  categorical.values <- rep(
    c("present", "other", "NA", NA_character_), length.out = ncol(input$object)
  )

  for (metadata.values in list(
    categorical.values, factor(categorical.values), rep(NA_character_, ncol(input$object))
  )) {
    input$object@meta.data$metadata_with_na <- metadata.values
    for (show.labels in c(FALSE, TRUE)) {
      input$show.labels <- show.labels
      output <- do.call(plotMetadata, input)
      plot <- output$plots[[1]]
      built.plot <- ggplot2::ggplot_build(plot)
      missing.cells <- is.na(plot$data$clusid)
      colour.scale <- built.plot$plot$scales$get_scales("colour")

      expect_equal(nrow(built.plot$data[[1]]), ncol(input$object))
      expect_identical(plot$data$clusid, as.character(metadata.values))
      expect_true(all(built.plot$data[[1]]$colour[missing.cells] == "grey50"))
      expect_true(anyNA(colour.scale$get_breaks()))
      if (any(!missing.cells)) {
        expect_false(any(built.plot$data[[1]]$colour[!missing.cells] == "grey50"))
      }
      expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
      expect_identical(output$object@meta.data, input$object@meta.data)
    }
  }
})

test_that("plotMetadata retains numeric NA cells in grey without losing the continuous scale", {
  input <- getParamPM("Chariou")
  input$metadata.to.plot <- "metadata_with_na"
  numeric.values <- as.numeric(seq_len(ncol(input$object)))
  numeric.values[seq.int(1L, length(numeric.values), by = 4L)] <- NA_real_

  for (metadata.values in list(numeric.values, rep(NA_real_, ncol(input$object)))) {
    input$object@meta.data$metadata_with_na <- metadata.values
    output <- do.call(plotMetadata, input)
    plot <- output$plots[[1]]
    built.plot <- ggplot2::ggplot_build(plot)
    missing.cells <- is.na(plot$data$clusid)
    colour.scale <- built.plot$plot$scales$get_scales("colour")
    missing.legend <- built.plot$plot$scales$get_scales("shape")

    expect_equal(nrow(built.plot$data[[1]]), ncol(input$object))
    expect_equal(sum(missing.cells), sum(is.na(metadata.values)))
    expect_true(all(built.plot$data[[1]]$colour[missing.cells] == "grey50"))
    expect_false(colour.scale$is_discrete())
    expect_equal(as.character(missing.legend$get_breaks()), "NA")
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
    expect_identical(output$object@meta.data, input$object@meta.data)
  }
})

test_that("plotMetadata applies dot.size to categorical and numeric metadata including NA cells", {
  input <- getParamPM("Chariou")
  input$object@meta.data$numeric_for_size <- as.numeric(seq_len(ncol(input$object)))
  input$object@meta.data$numeric_for_size[1] <- NA_real_
  input$metadata.to.plot <- c("orig.ident", "numeric_for_size")

  for (dot.size in c(0.2, 2.5)) {
    input$dot.size <- dot.size
    output <- do.call(plotMetadata, input)

    expect_length(output$plots, 2)
    for (plot in output$plots) {
      point.data <- ggplot2::ggplot_build(plot)$data[[1]]
      expect_equal(nrow(point.data), ncol(input$object))
      expect_true(all(point.data$size == dot.size))
    }
  }
})

test_that("plotMetadata preserves exact metadata column names", {
  chariou.data <- getParamPM("Chariou")
  original.metadata.columns <- colnames(chariou.data$object@meta.data)
  chariou.data$metadata.to.plot <- 'c("SCT_snn_res.2.4")'

  captured.messages <- character()
  captured.output <- capture.output(
    output <- withCallingHandlers(
      do.call(plotMetadata, chariou.data),
      message = function(condition) {
        captured.messages <<- c(captured.messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      }
    )
  )

  expect_length(output$plots, 1)
  expect_identical(colnames(output$object@meta.data), original.metadata.columns)
  expect_true(any(grepl("  - SCT_snn_res.2.4", captured.messages, fixed = TRUE)))
  expect_false(any(grepl("  - SCT_snn_res_2_4", captured.messages, fixed = TRUE)))
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
  output.without.labels <- do.call(plotMetadata, chariou.data)
  layer.geoms.without.labels <- vapply(
    output.without.labels$plots[[1]]$layers,
    function(layer) class(layer$geom)[1],
    character(1)
  )
  expect_false("GeomLabelRepel" %in% layer.geoms.without.labels)

  chariou.data$show.labels <- TRUE
  output.with.labels <- do.call(plotMetadata, chariou.data)
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

  output <- do.call(plotMetadata, chariou.data)
  plot <- output$plots[[1]]

  expect_identical(plot$theme$legend.position, "bottom")
  expect_equal(plot$theme$legend.text$size, ggplot2::rel(1.5))
})

test_that("Test Plot Metadata using TEC (Mouse) dataset", {
  tec.data <- getParamPM("TEC")
  output <- do.call(plotMetadata, tec.data)

  ggsave("output/TEC_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet1.png")
  ggsave("output/TEC_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet2.png")
  ggsave("output/TEC_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet3.png")
  ggsave("output/TEC_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet4.png")
  ggsave("output/TEC_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet5.png")
  ggsave("output/TEC_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet6.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using Chariou (Mouse) dataset", {
  chariou.data <- getParamPM("Chariou")
  output <- do.call(plotMetadata, chariou.data)

  ggsave("output/Chariou_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output/Chariou_plotmet.png", "Chariou_plotmet1.png")
  ggsave("output/Chariou_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output/Chariou_plotmet.png", "Chariou_plotmet2.png")
  ggsave("output/Chariou_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output/Chariou_plotmet.png", "Chariou_plotmet3.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using BRCA (Human) dataset", {
  brca.data <- getParamPM("BRCA")
  output <- do.call(plotMetadata, brca.data)

  ggsave("output/BRCA_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output/BRCA_plotmet.png", "BRCA_plotmet1.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output/BRCA_plotmet.png", "BRCA_plotmet2.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output/BRCA_plotmet.png", "BRCA_plotmet3.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output/BRCA_plotmet.png", "BRCA_plotmet4.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output/BRCA_plotmet.png", "BRCA_plotmet5.png")
  ggsave("output/BRCA_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output/BRCA_plotmet.png", "BRCA_plotmet6.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using NSCLCmulti (Human) dataset", {
  nsclc.multi.data <- getParamPM("nsclc-multi")
  output <- do.call(plotMetadata, nsclc.multi.data)

  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output/NSCLCmulti_plotmet.png", "NSCLCmulti_plotmet1.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output/NSCLCmulti_plotmet.png", "NSCLCmulti_plotmet2.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output/NSCLCmulti_plotmet.png", "NSCLCmulti_plotmet3.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output/NSCLCmulti_plotmet.png", "NSCLCmulti_plotmet4.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output/NSCLCmulti_plotmet.png", "NSCLCmulti_plotmet5.png")
  ggsave("output/NSCLCmulti_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output/NSCLCmulti_plotmet.png", "NSCLCmulti_plotmet6.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using PBMCsingle (Human) dataset", {
  pbmc.single.data <- getParamPM("pbmc-single")
  output <- do.call(plotMetadata, pbmc.single.data)

  ggsave("output/PBMCsingle_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output/PBMCsingle_plotmet.png", "PBMCsingle_plotmet1.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output/PBMCsingle_plotmet.png", "PBMCsingle_plotmet2.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output/PBMCsingle_plotmet.png", "PBMCsingle_plotmet3.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output/PBMCsingle_plotmet.png", "PBMCsingle_plotmet4.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output/PBMCsingle_plotmet.png", "PBMCsingle_plotmet5.png")
  ggsave("output/PBMCsingle_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output/PBMCsingle_plotmet.png", "PBMCsingle_plotmet6.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})

test_that("Test Plot Metadata using TEC (Mouse) dataset; UMAP", {
  tec.data <- getParamPM("TEC")
  tec.data$reduction.type <- "umap"
  output <- do.call(plotMetadata, tec.data)

  ggsave("output/TEC_plotmet.png", output$plot[[1]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet1.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[2]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet2.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[3]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet3.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[4]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet4.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[5]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet5.umap.png")
  ggsave("output/TEC_plotmet.png", output$plot[[6]], width = 10, height = 10)
  expect_snapshot_file("output/TEC_plotmet.png", "TEC_plotmet6.umap.png")

  expect_type(output, "list")
  expected.elements = c("object", "plots")
  expect_setequal(names(output), expected.elements)
})
