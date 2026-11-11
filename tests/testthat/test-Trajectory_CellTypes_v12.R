skip_if_packages_missing <- function(pkgs) {
  for (pkg in pkgs) {
    testthat::skip_if_not_installed(pkg)
  }
}

build_test_seurat <- function() {
  set.seed(42)
  counts <- matrix(
    rpois(200, lambda = 5),
    nrow = 20,
    dimnames = list(paste0("g", seq_len(20)), paste0("c", seq_len(10)))
  )
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts)
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- suppressWarnings(
    Seurat::FindVariableFeatures(
      seurat_obj,
      selection.method = "vst",
      nfeatures = 10,
      verbose = FALSE
    )
  )
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE)
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = 5, verbose = FALSE)
  seurat_obj <- Seurat::RunTSNE(
    seurat_obj,
    dims = 1:5,
    perplexity = 2,
    check_duplicates = FALSE,
    verbose = FALSE
  )
  seurat_obj <- Seurat::RunUMAP(
    seurat_obj,
    dims = 1:5,
    n.neighbors = 5,
    verbose = FALSE
  )
  clusters <- sample(letters[1:2], ncol(seurat_obj), replace = TRUE)
  names(clusters) <- colnames(seurat_obj)
  Seurat::AddMetaData(
    object = seurat_obj,
    metadata = data.frame(cluster = clusters, row.names = names(clusters))
  )
}

required_pkgs <- c(
  "Seurat", "scater", "scran", "scuttle", "TrajectoryUtils", "TSCAN",
  "ggplot2", "gridExtra", "mclust",
  "SingleCellExperiment", "SummarizedExperiment"
)

testthat::test_that("load_seurat_object handles objects, files and NIDAP inputs", {
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("g1", "g2"), c("c1", "c2"))
    )
  )
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(seurat_obj, rds_path)
  requested_file <- NULL
  input <- list(fileSystem = function() {
    list(get_path = function(name, type) {
      requested_file <<- list(name = name, type = type)
      rds_path
    })
  })
  loader <- load_seurat_object
  environment(loader) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = environment(loader)
  )

  testthat::expect_identical(loader(seurat_obj), seurat_obj)
  testthat::expect_identical(loader(rds_path), seurat_obj)
  testthat::expect_identical(loader(input), seurat_obj)
  testthat::expect_identical(
    requested_file,
    list(name = "seurat_object.rds", type = "r")
  )
})

testthat::test_that("load_seurat_object preserves errors for invalid inputs", {
  invalid_filesystem <- list(fileSystem = function() stop("Filesystem unavailable"))
  missing_file <- list(fileSystem = function() {
    list(get_path = function(name, type) tempfile(fileext = ".rds"))
  })
  expected_error <- "Provide a Seurat object or a path to an RDS file containing one."

  testthat::expect_error(load_seurat_object(NULL), expected_error, fixed = TRUE)
  testthat::expect_error(
    load_seurat_object(invalid_filesystem), expected_error, fixed = TRUE
  )
  testthat::expect_error(load_seurat_object(missing_file), expected_error, fixed = TRUE)
})


 testthat::test_that("Trajectory_CellTypes errors when clustering column missing", {
  skip_if_packages_missing(required_pkgs)
  seurat_obj <- build_test_seurat()
  metadata <- data.frame(other = rep("a", ncol(seurat_obj)), row.names = colnames(seurat_obj))
  testthat::expect_error(
    Trajectory_CellTypes(
      Seurat_Object = seurat_obj,
      MetaData = metadata,
      Clusters_to_Use = "cluster"
    ),
    "Column 'cluster' not found"
  )
})

 testthat::test_that("Trajectory_CellTypes validates argument types", {
  skip_if_packages_missing(required_pkgs)
  seurat_obj <- build_test_seurat()
  metadata <- data.frame(cluster = seurat_obj$cluster, row.names = colnames(seurat_obj))
  testthat::expect_error(
    Trajectory_CellTypes(
      Seurat_Object = seurat_obj,
      MetaData = metadata,
      Clusters_to_Use = 123
    ),
    "non-missing character scalar"
  )
})

 testthat::test_that("Trajectory_CellTypes runs on minimal input", {
  skip_if_packages_missing(required_pkgs)
  seurat_obj <- build_test_seurat()
  metadata <- data.frame(cluster = seurat_obj$cluster, row.names = colnames(seurat_obj))
  set.seed(7)
  testthat::expect_invisible(
    Trajectory_CellTypes(
      Seurat_Object = seurat_obj,
      MetaData = metadata,
      Clusters_to_Use = "cluster",
      Custom_Gene_List = "g1,g2"
    )
  )
})
