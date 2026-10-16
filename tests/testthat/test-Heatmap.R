test_that("Produce heatmap and return plot and filtered dataframe: TEC data",
{
  cr.object <- getParamHM("TEC")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)
  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                                 "TEC_heatmap.png")
})

test_that("Heatmap scaled vs unscaled", {
  cr.object <- getParamHM("TEC")
  output <- do.call(heatmapSC, cr.object)
  cr.object$scale.data <- FALSE
  output2 <- do.call(heatmapSC, cr.object)

  #compare scaled (a) vs. nonscaled data (b) to be different
  a <- rowMeans(as.data.frame.matrix(output$data)[, -1])
  b <- rowMeans(as.data.frame.matrix(output2$data)[, -1])

  expect_false(isTRUE(all.equal(a, b)))

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output2$plots),
                       "TEC_heatmap_unscaled.png")

})

test_that("Heatmap uses orig.ident without renaming metadata columns", {
  cr.object <- getParamHM("Chariou")
  metadata.names <- colnames(cr.object$object@meta.data)
  cr.object$metadata <- "orig.ident"

  output <- do.call(heatmapSC, cr.object)

  expect_setequal(names(output), c("plots", "data"))
  expect_identical(colnames(cr.object$object@meta.data), metadata.names)
  expect_true("orig.ident" %in% colnames(cr.object$object@meta.data))
  expect_false("orig_ident" %in% colnames(cr.object$object@meta.data))
})

test_that("Heatmap supports a custom sample name column", {
  cr.object <- getParamHM("Chariou")
  cr.object$object@meta.data[["sample.name_1.2"]] <-
    cr.object$object@meta.data[["orig.ident"]]
  metadata.names <- colnames(cr.object$object@meta.data)
  cr.object$sample.name.column <- "sample.name_1.2"
  selected.sample <- unique(cr.object$object@meta.data[["sample.name_1.2"]])[[1]]
  cr.object$sample.names <- selected.sample
  expected.cells <- rownames(cr.object$object@meta.data)[
    cr.object$object@meta.data[["sample.name_1.2"]] == selected.sample
  ]

  expect_message(
    output <- do.call(heatmapSC, cr.object),
    paste0(
      "Selected samples using metadata column 'sample.name_1.2': '",
      selected.sample, "'"
    ),
    fixed = TRUE
  )
  expect_message(
    do.call(heatmapSC, cr.object),
    paste0(
      "Displaying ", length(expected.cells), "/",
      nrow(cr.object$object@meta.data), " cells."
    ),
    fixed = TRUE
  )

  expect_setequal(names(output), c("plots", "data"))
  expect_setequal(colnames(output$data)[-1], expected.cells)
  expect_identical(colnames(cr.object$object@meta.data), metadata.names)
})

test_that("Heatmap uses all samples when sample selection is blank", {
  cr.object <- getParamHM("Chariou")
  cr.object$sample.names <- c("")

  expect_message(
    output <- do.call(heatmapSC, cr.object),
    "No samples specified. Using all samples.",
    fixed = TRUE
  )

  expect_setequal(colnames(output$data)[-1], rownames(cr.object$object@meta.data))
})

test_that("Heatmap preserves punctuation in sample names", {
  cr.object <- getParamHM("Chariou")
  sample.column <- "sample.name_1.2"
  sample.values <- rep(
    c("sample.1-A", "sample.2-B"),
    length.out = nrow(cr.object$object@meta.data)
  )
  cr.object$object@meta.data[[sample.column]] <- sample.values
  metadata.names <- colnames(cr.object$object@meta.data)
  cr.object$sample.name.column <- sample.column
  cr.object$sample.names <- "sample.1-A"
  expected.cells <- rownames(cr.object$object@meta.data)[
    cr.object$object@meta.data[[sample.column]] == "sample.1-A"
  ]

  output <- do.call(heatmapSC, cr.object)

  expect_setequal(colnames(output$data)[-1], expected.cells)
  expect_identical(colnames(cr.object$object@meta.data), metadata.names)
})

test_that("Heatmap defaults blank sample name columns to orig.ident", {
  cr.object <- getParamHM("Chariou")
  cr.object$sample.name.column <- ""

  output <- do.call(heatmapSC, cr.object)

  expect_setequal(names(output), c("plots", "data"))
})

test_that("Heatmap warns when the default sample name column is absent", {
  cr.object <- getParamHM("Chariou")
  cr.object$sample.name.column <- ""
  cr.object$object@meta.data[["orig.ident"]] <- NULL

  expect_warning(
    expect_error(
      do.call(heatmapSC, cr.object),
      "Could not find requested sample name column: orig.ident"
    ),
    "Default sample name column 'orig.ident' was not found in Seurat metadata."
  )
})

test_that("Heatmap rejects invalid sample name columns", {
  cr.object <- getParamHM("Chariou")
  cr.object$sample.name.column <- c("orig.ident", "Phase")

  expect_error(
    do.call(heatmapSC, cr.object),
    "sample.name.column must be a single, non-empty metadata column name."
  )

  cr.object$sample.name.column <- "not_a_metadata_column"
  expect_error(
    do.call(heatmapSC, cr.object),
    "Could not find requested sample name column: not_a_metadata_column"
  )
})

test_that("Heatmap identifies unavailable requested sample names", {
  cr.object <- getParamHM("Chariou")
  available.samples <- sort(unique(as.character(
    cr.object$object@meta.data[["orig.ident"]]
  )))
  selected.sample <- available.samples[[1]]
  unavailable.sample <- "not_a_sample"
  cr.object$sample.names <- c(selected.sample, unavailable.sample)
  expected.message <- paste0(
    "Requested sample name(s) not found in metadata column 'orig.ident': '",
    unavailable.sample,
    "'.\nAvailable sample names: ",
    paste(shQuote(available.samples), collapse = ", ")
  )

  expect_error(
    do.call(heatmapSC, cr.object),
    expected.message,
    fixed = TRUE
  )
})

test_that("Heatmap supports a single transcript", {
  cr.object <- getParamHM("Chariou")
  cr.object$transcripts <- cr.object$transcripts[[1]]
  cr.object$add.gene.or.protein <- FALSE
  cr.object$rna.annotations <- NULL

  output <- do.call(heatmapSC, cr.object)

  expect_identical(output$data$gene, cr.object$transcripts)
})

test_that("Heatmap supports a single protein", {
  cr.object <- getParamHM("nsclc-multi")
  cr.object$transcripts <- cr.object$transcripts[[1]]
  cr.object$proteins <- cr.object$proteins[[1]]
  cr.object$add.gene.or.protein <- FALSE
  cr.object$rna.annotations <- NULL
  cr.object$protein.annotations <- NULL

  output <- do.call(heatmapSC, cr.object)

  expect_setequal(
    output$data$gene,
    c(cr.object$transcripts, paste0(cr.object$proteins, "_Prot"))
  )
})

test_that("Heatmap run with bad gene name", {
  cr.object <- getParamHM("TEC")

  #Spike in with bad gene name
  cr.object$transcripts <- c("badgene", cr.object$transcripts)
  missing.transcript.message <- paste0(
    "There are 1 gene(s) absent from dataset: 'badgene'.\n",
    "          Possible reasons are that gene is not official gene symbol or gene ",
    "is not highly expressed and has been filtered."
  )

  expect_no_warning(
    logs <- capture.output(
      do.call(heatmapSC, cr.object),
      type = "message"
    )
  )
  log.text <- paste(logs, collapse = "\n")
  expect_true(grepl(missing.transcript.message, log.text, fixed = TRUE))
  expect_lt(
    grep("Transcripts found in", logs, fixed = TRUE)[[1]],
    grep("There are 1 gene(s) absent from dataset", logs, fixed = TRUE)[[1]]
  )
})

test_that("Heatmap reports missing proteins and retains available proteins", {
  cr.object <- getParamHM("nsclc-multi")
  cr.object$proteins <- c("CD3", "not_a_protein")
  cr.object$protein.annotations <- NULL

  expect_no_warning(
    logs <- capture.output(
      output <- do.call(heatmapSC, cr.object),
      type = "message"
    )
  )
  missing.protein.message <-
    "Proteins not found in Protein assay: not_a_protein. They will be excluded."
  expect_true(any(grepl(missing.protein.message, logs, fixed = TRUE)))
  expect_lt(
    grep("Proteins found in", logs, fixed = TRUE)[[1]],
    grep(missing.protein.message, logs, fixed = TRUE)[[1]]
  )

  expect_true("CD3_Prot" %in% output$data$gene)
  expect_false("not_a_protein_Prot" %in% output$data$gene)
})

test_that("Heatmap requires the selected transcript assay", {
  cr.object <- getParamHM("TEC")
  cr.object$use.assay <- "Harmony"

  expect_error(
    do.call(heatmapSC, cr.object),
    "Transcript assay 'Harmony' was not found in the Seurat object."
  )
})

test_that("Heatmap requires a CITE-seq Protein assay for proteins", {
  cr.object <- getParamHM("TEC")
  cr.object$transcripts <- character()
  cr.object$proteins <- "CD3"

  expect_error(
    do.call(heatmapSC, cr.object),
    "Protein assay was not found.",
    fixed = TRUE
  )
})

test_that("Heatmap run with duplicate gene names", {
  cr.object <- getParamHM("TEC")

  #Spike in with bad gene name
  cr.object$transcripts <- c("Map7", cr.object$transcripts)

  expect_warning(do.call(heatmapSC, cr.object),
                 "The following duplicate genes were removed: Map7")

})

test_that("Heatmap run with error for missing all genes", {
  cr.object <- getParamHM("TEC")
  cr.object$transcripts <-
    c("APOE", "ARG1", "CD38", "CD3D", "CD3E", "CD3G")
  missing.transcript.message <- paste0(
    "There are 6 gene(s) absent from dataset: 'APOE', 'ARG1', 'CD38', 'CD3D', ",
    "'CD3E', 'CD3G'.\n          Possible reasons are that gene is not official ",
    "gene symbol or gene is not highly expressed and has been filtered."
  )
  expect_no_warning(
    expect_message(
      expect_error(
        do.call(heatmapSC, cr.object),
        "No requested transcripts or proteins were found in the selected assays."
      ),
      missing.transcript.message,
      fixed = TRUE
    )
  )

})

test_that("Produce heatmap - Chariou data", {
  cr.object <- getParamHM("Chariou")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "Chariou_heatmap.png")
})

test_that("Chariou with no additional protein/transcript annotations", {
  cr.object <- getParamHM("Chariou2")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "Chariou_heatmap2.png")
})

test_that("Produce heatmap - Chariou with no transcripts/proteins", {
  cr.object <- getParamHM("Chariou3")

  expect_error(do.call(heatmapSC, cr.object),
               "At least 1 transcript and/or protein is needed for plotting")
})

test_that("Produce heatmap - PBMC single data", {
  cr.object <- getParamHM("pbmc-single")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "pbmc-single_heatmap.png")
})

test_that("Produce heatmap with protein using NSCLC multi data", {
  cr.object <- getParamHM("nsclc-multi")
  cr.object$legend.font = 3
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "nsclc-multi_heatmap.png")
})

test_that("Produce heatmap with filtered dataframe - BRCA data", {
  cr.object <- getParamHM("BRCA")
  output <- do.call(heatmapSC, cr.object)

  expect_type(output, "list")
  expected.elements = c("plots", "data")
  expect_setequal(names(output), expected.elements)

  skip_on_ci()
  expect_snapshot_file(.drawHeatPng(output$plots),
                       "BRCA_heatmap.png")
})